/// @file IndexedBVH.cpp
/// @brief Device-transfer support for the compact IndexedBVH representation.
///
/// Started June 2026 by Sandro Wenzel (sandro.wenzel@cern.ch) and Tristan Wenzel
/// (tristan.wenzel@cern.ch).

#include "VecGeom/base/IndexedBVH.h"

#ifdef VECGEOM_CUDA_INTERFACE
#include "VecGeom/backend/cuda/Interface.h"
#include <cassert>
#include <cstring>
#include <stdexcept>
#endif

namespace vecgeom {
namespace cuda {
template <typename Real_t>
void MaterializeIndexedBVHDeviceImage(void *arena, IndexedBVH<Real_t> *target);
template <typename Real_t>
void MaterializeIndexedBVHDeviceImages(void *arena, std::size_t count, std::size_t objectStride);
} // namespace cuda
inline namespace VECGEOM_IMPL_NAMESPACE {

#ifdef VECGEOM_CUDA_INTERFACE
namespace {

std::size_t AlignDeviceOffset(std::size_t value, std::size_t alignment)
{
  assert(alignment > 0 && (alignment & (alignment - 1)) == 0);
  return (value + alignment - 1) & ~(alignment - 1);
}

template <typename T>
std::size_t ReserveDeviceArray(std::size_t &nextOffset, std::size_t count)
{
  if (count == 0) return 0;
  nextOffset        = AlignDeviceOffset(nextOffset, alignof(T));
  const auto offset = nextOffset;
  nextOffset += count * sizeof(T);
  return offset;
}

template <typename T>
void PackDeviceArray(char *arena, std::size_t arenaSize, std::size_t offset, T const *source, std::size_t count)
{
  if (count == 0) return;
  assert(offset > 0 && source != nullptr && offset + count * sizeof(T) <= arenaSize);
  std::memcpy(arena + offset, source, count * sizeof(T));
}

} // namespace

template <typename Real_t>
typename IndexedBVH<Real_t>::DeviceImageLayout IndexedBVH<Real_t>::MakeDeviceImageLayout(std::size_t &nextOffset) const
{
  DeviceImageLayout layout;
  layout.rootId              = fRootId;
  layout.nPrim               = fNPrim;
  layout.nDaughter           = fNDaughter;
  layout.nNodes              = fNNodes;
  layout.nodesOffset         = ReserveDeviceArray<Node>(nextOffset, fNNodes);
  layout.primIdOffset        = ReserveDeviceArray<int>(nextOffset, fNPrim);
  layout.aabbsOffset         = ReserveDeviceArray<AABB<Real_t>>(nextOffset, fNPrim);
  layout.boxToDaughterOffset = fBoxToDaughter ? ReserveDeviceArray<int>(nextOffset, fNPrim) : 0;
  return layout;
}

template <typename Real_t>
void IndexedBVH<Real_t>::PackDeviceImage(char *arena, std::size_t arenaSize, std::size_t objectOffset,
                                         DeviceImageLayout const &layout) const
{
  assert(arena != nullptr && objectOffset + sizeof(layout) <= arenaSize);
  assert(layout.rootId == fRootId && layout.nPrim == fNPrim && layout.nDaughter == fNDaughter &&
         layout.nNodes == fNNodes);
  std::memcpy(arena + objectOffset, &layout, sizeof(layout));
  PackDeviceArray(arena, arenaSize, layout.nodesOffset, fNodes, fNNodes);
  PackDeviceArray(arena, arenaSize, layout.primIdOffset, fPrimId, fNPrim);
  PackDeviceArray(arena, arenaSize, layout.aabbsOffset, fAABBs, fNPrim);
  PackDeviceArray(arena, arenaSize, layout.boxToDaughterOffset, fBoxToDaughter, fBoxToDaughter ? fNPrim : 0);
}

template <typename Real_t>
DevicePtr<cuda::IndexedBVH<Real_t>> IndexedBVH<Real_t>::CopyToGpu(void *addr) const
{
  if (!addr) throw std::logic_error("Cannot copy IndexedBVH into a null pointer!");

  std::size_t arenaSize = sizeof(DeviceImageLayout);
  const auto layout     = MakeDeviceImageLayout(arenaSize);
  std::vector<char> hostArena(arenaSize, 0);
  PackDeviceImage(hostArena.data(), hostArena.size(), 0, layout);

  char *deviceArena = nullptr;
  VECGEOM_DEVICE_API_CALL(Malloc((void **)&deviceArena, arenaSize));
  VECGEOM_DEVICE_API_CALL(
      Memcpy(deviceArena, hostArena.data(), arenaSize, VECGEOM_DEVICE_API_SYMBOL(MemcpyHostToDevice)));
  auto *target = reinterpret_cast<cuda::IndexedBVH<Real_t> *>(addr);
  cuda::MaterializeIndexedBVHDeviceImage<Real_t>(deviceArena, target);
  return DevicePtr<cuda::IndexedBVH<Real_t>>(target);
}

template DevicePtr<cuda::IndexedBVH<float>> IndexedBVH<float>::CopyToGpu(void *) const;
template DevicePtr<cuda::IndexedBVH<double>> IndexedBVH<double>::CopyToGpu(void *) const;
template IndexedBVH<float>::DeviceImageLayout IndexedBVH<float>::MakeDeviceImageLayout(std::size_t &) const;
template IndexedBVH<double>::DeviceImageLayout IndexedBVH<double>::MakeDeviceImageLayout(std::size_t &) const;
template void IndexedBVH<float>::PackDeviceImage(char *, std::size_t, std::size_t, DeviceImageLayout const &) const;
template void IndexedBVH<double>::PackDeviceImage(char *, std::size_t, std::size_t, DeviceImageLayout const &) const;
#endif // VECGEOM_CUDA_INTERFACE

} // namespace VECGEOM_IMPL_NAMESPACE

#ifdef VECGEOM_CUDA_INTERFACE
template <typename Real_t>
cuda::IndexedBVH<Real_t> *CopyManyIndexedBVHToGpu(void const *rawHostBVHs, std::size_t count)
{
  using HostBVH = cxx::IndexedBVH<Real_t>;
  if (count == 0) return nullptr;
  if (!rawHostBVHs) throw std::logic_error("Cannot copy a null IndexedBVH array to the GPU!");

  auto const *hostBVHs = static_cast<HostBVH const *>(rawHostBVHs);
  std::vector<typename HostBVH::DeviceImageLayout> layouts(count);
  std::size_t arenaSize = count * sizeof(HostBVH);
  for (std::size_t index = 0; index < count; ++index)
    layouts[index] = hostBVHs[index].MakeDeviceImageLayout(arenaSize);

  std::vector<char> hostArena(arenaSize, 0);
  for (std::size_t index = 0; index < count; ++index)
    hostBVHs[index].PackDeviceImage(hostArena.data(), hostArena.size(), index * sizeof(HostBVH), layouts[index]);

  char *deviceArena = nullptr;
  VECGEOM_DEVICE_API_CALL(Malloc((void **)&deviceArena, arenaSize));
  VECGEOM_DEVICE_API_CALL(
      Memcpy(deviceArena, hostArena.data(), arenaSize, VECGEOM_DEVICE_API_SYMBOL(MemcpyHostToDevice)));
  cuda::MaterializeIndexedBVHDeviceImages<Real_t>(deviceArena, count, sizeof(HostBVH));
  return reinterpret_cast<cuda::IndexedBVH<Real_t> *>(deviceArena);
}

template cuda::IndexedBVH<float> *CopyManyIndexedBVHToGpu<float>(void const *, std::size_t);
template cuda::IndexedBVH<double> *CopyManyIndexedBVHToGpu<double>(void const *, std::size_t);
#endif

#ifdef VECCORE_CUDA
namespace cuda {

template <typename Real_t>
__global__ void MaterializeIndexedBVHDeviceImageKernel(char *arena, IndexedBVH<Real_t> *target)
{
  IndexedBVH<Real_t>::MaterializeDeviceImage(arena, 0, target);
}

template <typename Real_t>
void MaterializeIndexedBVHDeviceImage(void *arena, IndexedBVH<Real_t> *target)
{
  MaterializeIndexedBVHDeviceImageKernel<Real_t><<<1, 1>>>(reinterpret_cast<char *>(arena), target);
  VECGEOM_DEVICE_API_CALL(GetLastError());
  VECGEOM_DEVICE_API_CALL(DeviceSynchronize());
}

template <typename Real_t>
__global__ void MaterializeIndexedBVHDeviceImagesKernel(char *arena, std::size_t count, std::size_t objectStride)
{
  for (std::size_t index = blockIdx.x * blockDim.x + threadIdx.x; index < count; index += blockDim.x * gridDim.x)
    IndexedBVH<Real_t>::MaterializeDeviceImage(arena, index * objectStride);
}

template <typename Real_t>
void MaterializeIndexedBVHDeviceImages(void *arena, std::size_t count, std::size_t objectStride)
{
  VECGEOM_ASSERT(objectStride == sizeof(IndexedBVH<Real_t>));
  constexpr int blockSize = 256;
  const int gridSize      = static_cast<int>((count + blockSize - 1) / blockSize);
  MaterializeIndexedBVHDeviceImagesKernel<Real_t>
      <<<gridSize, blockSize>>>(reinterpret_cast<char *>(arena), count, objectStride);
  VECGEOM_DEVICE_API_CALL(GetLastError());
  VECGEOM_DEVICE_API_CALL(DeviceSynchronize());
}

template void MaterializeIndexedBVHDeviceImage<float>(void *, IndexedBVH<float> *);
template void MaterializeIndexedBVHDeviceImage<double>(void *, IndexedBVH<double> *);
template void MaterializeIndexedBVHDeviceImages<float>(void *, std::size_t, std::size_t);
template void MaterializeIndexedBVHDeviceImages<double>(void *, std::size_t, std::size_t);

} // namespace cuda
#endif

} // namespace vecgeom
