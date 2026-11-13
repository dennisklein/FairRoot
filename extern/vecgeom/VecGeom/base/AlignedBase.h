/// \file AlignedBase.h
/// \author Johannes de Fine Licht (johannes.definelicht@cern.ch)

#ifndef VECGEOM_BASE_ALIGNEDBASE_H_
#define VECGEOM_BASE_ALIGNEDBASE_H_

#include "VecGeom/base/Global.h"

#include <new>

#if VECGEOM_ALIGNEDBASE_ALIGNMENT > 1 && defined(__cpp_aligned_new)
// C++17 routes placement-new for over-aligned types through the aligned
// placement form. Provide the non-allocating overload so raw-buffer
// construction keeps working for host and CUDA helper layouts.
VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE void *operator new(std::size_t, std::align_val_t, void *p) noexcept
{
  return p;
}

VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE void operator delete(void *, std::align_val_t, void *) noexcept {}
#endif

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE(class AlignedBase;);
VECGEOM_DEVICE_DECLARE_CONV(class, AlignedBase);

inline namespace VECGEOM_IMPL_NAMESPACE {

#if VECGEOM_ALIGNEDBASE_ALIGNMENT > 1
constexpr size_t kAlignedBaseAlignment = VECGEOM_ALIGNEDBASE_ALIGNMENT;

#ifndef VECCORE_CUDA
static_assert(kAlignedBaseAlignment == alignof(vecgeom::VectorBackend::Real_v),
              "Configured AlignedBase alignment must match the host vector backend.");
#endif

// Keep host/device helper layout identical while preserving Vc-friendly stride
// alignment for raw-copied arrays such as Vector3D pools.
class alignas(kAlignedBaseAlignment) AlignedBase {
public:
#ifndef VECCORE_CUDA
  VECGEOM_FORCE_INLINE
  void *operator new(size_t size) { return vecCore::AlignedAlloc(kAlignedBaseAlignment, size); }

  VECGEOM_FORCE_INLINE
  void *operator new(size_t, void *p) { return p; }

  VECGEOM_FORCE_INLINE
  void *operator new[](size_t size) { return vecCore::AlignedAlloc(kAlignedBaseAlignment, size); }

  VECGEOM_FORCE_INLINE
  void *operator new[](size_t, void *p) { return p; }

  VECGEOM_FORCE_INLINE
  void operator delete(void *ptr, size_t) { vecCore::AlignedFree(ptr); }

  VECGEOM_FORCE_INLINE
  void operator delete(void *, void *) {}

  VECGEOM_FORCE_INLINE
  void operator delete[](void *ptr, size_t) { vecCore::AlignedFree(ptr); }

  VECGEOM_FORCE_INLINE
  void operator delete[](void *, void *) {}
#endif
};
#elif !defined(VECCORE_CUDA)
class AlignedBase {

public:
  VECGEOM_FORCE_INLINE
  void *operator new(size_t size) { return vecCore::AlignedAlloc(kAlignmentBoundary, size); }

  VECGEOM_FORCE_INLINE
  void *operator new(size_t, void *p) { return p; }

  VECGEOM_FORCE_INLINE
  void *operator new[](size_t size) { return vecCore::AlignedAlloc(kAlignmentBoundary, size); }

  VECGEOM_FORCE_INLINE
  void *operator new[](size_t, void *p) { return p; }

  VECGEOM_FORCE_INLINE
  void operator delete(void *ptr, size_t) { vecCore::AlignedFree(ptr); }

  VECGEOM_FORCE_INLINE
  void operator delete(void *, void *) {}

  VECGEOM_FORCE_INLINE
  void operator delete[](void *ptr, size_t) { vecCore::AlignedFree(ptr); }

  VECGEOM_FORCE_INLINE
  void operator delete[](void *, void *) {}
};
#else
class AlignedBase {};
#endif
} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif // VECGEOM_BASE_ALIGNEDBASE_H_
