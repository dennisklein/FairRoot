/// \file BVHManager.cu
/// \author Guilherme Amadio

#include "VecGeom/management/BVHManager.h"
#include "VecGeom/backend/cuda/Interface.h"
#include "VecGeom/navigation/BVHSafetyEstimator.h"

#include "VecGeom/base/Assert.h"

namespace vecgeom {
inline namespace cuda {
template <typename Real_t>
__global__ void InitializeDeviceIndexedBVHs(char *arena, size_t n)
{
  using DeviceBVH = IndexedBVH<Real_t>;
  using Layout    = typename DeviceBVH::DeviceImageLayout;
  static_assert(sizeof(Layout) <= sizeof(DeviceBVH), "IndexedBVH device layout does not fit in an object slot");
  static_assert(alignof(Layout) <= alignof(DeviceBVH), "IndexedBVH device layout requires stronger alignment");

  for (size_t id = blockIdx.x * blockDim.x + threadIdx.x; id < n; id += blockDim.x * gridDim.x)
    DeviceBVH::MaterializeDeviceImage(arena, id * sizeof(DeviceBVH));
}

/** Construct and publish IndexedBVH views over one packed device arena. */
template <typename Real_t>
inline IndexedBVH<Real_t> *InitializeDeviceIndexedBVHBuffer(void *arena, size_t n, size_t objectStride,
                                                            size_t layoutSize)
{
  using Layout = typename IndexedBVH<Real_t>::DeviceImageLayout;
  VECGEOM_ASSERT(arena != nullptr && n > 0);
  VECGEOM_ASSERT(objectStride == sizeof(IndexedBVH<Real_t>));
  VECGEOM_ASSERT(layoutSize == sizeof(Layout));
  auto *ptr               = reinterpret_cast<IndexedBVH<Real_t> *>(arena);
  constexpr int blockSize = 256;
  const int gridSize      = static_cast<int>((n + blockSize - 1) / blockSize);
  InitializeDeviceIndexedBVHs<Real_t><<<gridSize, blockSize>>>(reinterpret_cast<char *>(arena), n);
  VECGEOM_DEVICE_API_CALL(GetLastError());
  VECGEOM_DEVICE_API_CALL(MemcpyToSymbol(dBVH<Real_t>, &ptr, sizeof(ptr)));
  VECGEOM_DEVICE_API_CALL(DeviceSynchronize());
  return ptr;
}

/**
 * @brief Read back the published IndexedBVH device-array address.
 * @tparam Real_t Hierarchy storage precision.
 * @return Device-array base address stored in `dBVH<Real_t>`.
 */
template <typename Real_t>
inline IndexedBVH<Real_t> *GetDeviceIndexedBVHBuffer()
{
  IndexedBVH<Real_t> *ptr = nullptr;
  VECGEOM_DEVICE_API_CALL(MemcpyFromSymbol(&ptr, dBVH<Real_t>, sizeof(ptr)));
  return ptr;
}

// Temporary hack (used already in LogicalVolume.cpp) implementing the Instance functionality
// on device for BVHSafetyEstimator and BVHNavigatorV in the absence of the corresponding
// implementation files
VECCORE_ATT_DEVICE
BVHSafetyEstimator *gBVHSafetyEstimator = nullptr;

VECCORE_ATT_DEVICE
VSafetyEstimator *BVHSafetyEstimator::Instance()
{
  if (gBVHSafetyEstimator == nullptr) gBVHSafetyEstimator = new BVHSafetyEstimator();
  return gBVHSafetyEstimator;
}

template IndexedBVH<float> *InitializeDeviceIndexedBVHBuffer<float>(void *, size_t, size_t, size_t);
template IndexedBVH<double> *InitializeDeviceIndexedBVHBuffer<double>(void *, size_t, size_t, size_t);

template IndexedBVH<float> *GetDeviceIndexedBVHBuffer<float>();
template IndexedBVH<double> *GetDeviceIndexedBVHBuffer<double>();

} // namespace cuda
} // namespace vecgeom
