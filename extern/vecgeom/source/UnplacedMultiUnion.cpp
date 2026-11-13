#include "VecGeom/volumes/UnplacedMultiUnion.h"
#include "VecGeom/volumes/SpecializedMultiUnion.h"
#include "VecGeom/base/RNG.h"
#include <stdio.h>

#ifdef VECGEOM_CUDA_INTERFACE
#include "VecGeom/backend/cuda/Interface.h"
#include "VecGeom/management/CudaManager.h"

#include <vector>
#endif

namespace vecgeom {

#ifdef VECGEOM_CUDA_INTERFACE
namespace cuda {
template <typename Real_t>
void MaterializeIndexedBVHDeviceImages(void *arena, std::size_t count, std::size_t objectStride);
} // namespace cuda
#endif

inline namespace VECGEOM_IMPL_NAMESPACE {

Precision UnplacedMultiUnion::Capacity() const
{
  if (fMultiUnion.fCapacity >= 0.) {
    return fMultiUnion.fCapacity;
  }
  // Sample in the solid extent and estimate capacity by counting how many points are
  // sampled inside
  const size_t nsamples = 100000;
  fMultiUnion.fCapacity = EstimateCapacity(nsamples);
  return fMultiUnion.fCapacity;
}

Precision UnplacedMultiUnion::SurfaceArea() const
{
  if (fMultiUnion.fSurfaceArea >= 0.) {
    return fMultiUnion.fSurfaceArea;
  }
  // Sample points on components, count how many are on the solid surface to correct
  // the surface area of each component
  const size_t nsamples    = 10000;
  fMultiUnion.fSurfaceArea = 0.;
  for (size_t i = 0; i < GetNumberOfSolids(); ++i) {
    size_t nsurf = 0;
    for (size_t ip = 0; ip < nsamples; ++ip) {
      Vector3D<Precision> point =
          GetNode(i)->GetTransformation()->InverseTransform(GetNode(i)->GetUnplacedVolume()->SamplePointOnSurface());
      if (Inside(point) == vecgeom::kSurface) nsurf++;
    }
    fMultiUnion.fSurfaceArea += GetNode(i)->SurfaceArea() * nsurf / nsamples;
  }
  return fMultiUnion.fSurfaceArea;
}

Vector3D<Precision> UnplacedMultiUnion::SamplePointOnSurface() const
{
  // Select a random component
  auto counter                = 0;
  VPlacedVolume const *volume = nullptr;
  Vector3D<Precision> point;
  size_t id = 0;
  do {
    if (counter == 0) {
      id     = (size_t)RNG::Instance().uniform(0., GetNumberOfSolids());
      volume = GetNode(id);
    }
    point   = volume->GetTransformation()->InverseTransform(volume->GetUnplacedVolume()->SamplePointOnSurface());
    counter = (counter + 1) % 1000;
  } while (Inside(point) != vecgeom::kSurface);
  return point;
}

bool UnplacedMultiUnion::Normal(Vector3D<Precision> const &point, Vector3D<Precision> &normal) const
{
  // Compute normal to solid in a point
  bool valid = false;
  normal     = MultiUnionImplementation::NormalKernel<Precision>(fMultiUnion, point, valid);
  return valid;
}

#ifndef VECCORE_CUDA
VPlacedVolume *UnplacedMultiUnion::Create(LogicalVolume const *const logical_volume,
                                          Transformation3D const *const transformation, VPlacedVolume *const placement)
{
  if (placement) {
    new (placement) SpecializedMultiUnion(logical_volume, transformation);
    return placement;
  }
  return new SpecializedMultiUnion(logical_volume, transformation);
}

VPlacedVolume *UnplacedMultiUnion::SpecializedVolume(LogicalVolume const *const volume,
                                                     Transformation3D const *const transformation,
                                                     VPlacedVolume *const placement) const
{
  return VolumeFactory::CreateByTransformation<UnplacedMultiUnion>(volume, transformation, placement);
}
#else

VECCORE_ATT_DEVICE VPlacedVolume *UnplacedMultiUnion::Create(LogicalVolume const *const logical_volume,
                                                             Transformation3D const *const transformation, const int id,
                                                             const int copy_no, const int child_id,
                                                             VPlacedVolume *const placement)
{
  if (placement) {
    new (placement) SpecializedMultiUnion(logical_volume, transformation, id, copy_no, child_id);
    return placement;
  }
  return new SpecializedMultiUnion(logical_volume, transformation, id, copy_no, child_id);
}

VECCORE_ATT_DEVICE
VPlacedVolume *UnplacedMultiUnion::SpecializedVolume(LogicalVolume const *const volume,
                                                     Transformation3D const *const transformation, const int id,
                                                     const int copy_no, const int child_id,
                                                     VPlacedVolume *const placement) const
{
  return VolumeFactory::CreateByTransformation<UnplacedMultiUnion>(volume, transformation, id, copy_no, child_id,
                                                                   placement);
}

#endif

#ifdef VECGEOM_CUDA_INTERFACE

DevicePtr<cuda::VUnplacedVolume> UnplacedMultiUnion::CopyToGpu(DevicePtr<cuda::VUnplacedVolume> const in_gpu_ptr) const
{
  const size_t nvolumes = GetNumberOfSolids();
  VECGEOM_VALIDATE(nvolumes > 0 && fMultiUnion.fNavHelper && fMultiUnion.fNeighbourOffsets,
                   << "Cannot copy an open or empty MultiUnion to the GPU!");

  auto &manager = CudaManager::Instance();

  std::vector<cuda::VPlacedVolume const *> componentPointers(nvolumes);
  for (size_t i = 0; i < nvolumes; ++i)
    componentPointers[i] = manager.LookupPlaced(GetNode(i)).GetPtr();
  auto componentStorage = manager.AllocateDeviceBuffer(nvolumes * sizeof(componentPointers[0]));
  vecgeom::CopyToGpu(reinterpret_cast<char const *>(componentPointers.data()), componentStorage.GetPtr(),
                     nvolumes * sizeof(componentPointers[0]));
  auto **gpuComponents = reinterpret_cast<cuda::VPlacedVolume const **>(componentStorage.GetPtr());

  const size_t noffsets = nvolumes + 1;
  auto offsetStorage    = manager.AllocateDeviceBuffer(noffsets * sizeof(size_t));
  vecgeom::CopyToGpu(reinterpret_cast<char const *>(fMultiUnion.fNeighbourOffsets), offsetStorage.GetPtr(),
                     noffsets * sizeof(size_t));
  auto *gpuOffsets = reinterpret_cast<size_t *>(offsetStorage.GetPtr());

  const size_t nneighbours = fMultiUnion.fNeighbourOffsets[nvolumes];
  size_t *gpuNeighbours    = nullptr;
  if (nneighbours > 0) {
    auto neighbourStorage = manager.AllocateDeviceBuffer(nneighbours * sizeof(size_t));
    vecgeom::CopyToGpu(reinterpret_cast<char const *>(fMultiUnion.fNeighbours), neighbourStorage.GetPtr(),
                       nneighbours * sizeof(size_t));
    gpuNeighbours = reinterpret_cast<size_t *>(neighbourStorage.GetPtr());
  }

  using HostBVH            = MultiUnionStruct::BVHStructure;
  std::size_t bvhArenaSize = sizeof(HostBVH);
  const auto bvhLayout     = fMultiUnion.fNavHelper->MakeDeviceImageLayout(bvhArenaSize);
  std::vector<char> hostBVHArena(bvhArenaSize, 0);
  fMultiUnion.fNavHelper->PackDeviceImage(hostBVHArena.data(), hostBVHArena.size(), 0, bvhLayout);
  auto bvhStorage = manager.AllocateDeviceBuffer(bvhArenaSize);
  vecgeom::CopyToGpu(hostBVHArena.data(), bvhStorage.GetPtr(), bvhArenaSize);
  cuda::MaterializeIndexedBVHDeviceImages<MultiUnionStruct::BVHReal_t>(bvhStorage.GetPtr(), 1, sizeof(HostBVH));
  auto *gpuBVH = reinterpret_cast<cuda::IndexedBVH<MultiUnionStruct::BVHReal_t> *>(bvhStorage.GetPtr());

  return CopyToGpuImpl<UnplacedMultiUnion>(in_gpu_ptr, nvolumes, gpuComponents, gpuBVH, gpuOffsets, gpuNeighbours,
                                           fMultiUnion.fMinExtent.x(), fMultiUnion.fMinExtent.y(),
                                           fMultiUnion.fMinExtent.z(), fMultiUnion.fMaxExtent.x(),
                                           fMultiUnion.fMaxExtent.y(), fMultiUnion.fMaxExtent.z());
}

DevicePtr<cuda::VUnplacedVolume> UnplacedMultiUnion::CopyToGpu() const { return CopyToGpuImpl<UnplacedMultiUnion>(); }

#endif // VECGEOM_CUDA_INTERFACE

} // namespace VECGEOM_IMPL_NAMESPACE

#ifdef VECCORE_CUDA

namespace cxx {

template size_t DevicePtr<cuda::UnplacedMultiUnion>::SizeOf();
template void DevicePtr<cuda::UnplacedMultiUnion>::Construct() const;
template void DevicePtr<cuda::UnplacedMultiUnion>::Construct(size_t, cuda::VPlacedVolume const **,
                                                             cuda::IndexedBVH<float> *, size_t *, size_t *, Precision,
                                                             Precision, Precision, Precision, Precision,
                                                             Precision) const;

} // namespace cxx
#endif

} // namespace vecgeom
