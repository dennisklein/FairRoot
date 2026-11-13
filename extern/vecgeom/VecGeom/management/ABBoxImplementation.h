//===-- management/ABBoxImplementation.h -------------------------*- C++ -*-===//
//===----------------------------------------------------------------------===//
/// @file ABBoxImplementation.h
/// @brief Aligned bounding-box helper kernels for navigation acceleration.

#ifndef VECGEOM_MANAGEMENT_ABBOXIMPLEMENTATION_H_
#define VECGEOM_MANAGEMENT_ABBOXIMPLEMENTATION_H_

#include "VecGeom/base/Vector3D.h"
#include "VecGeom/volumes/kernel/GenericKernels.h"
#include <VecCore/VecCore>

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

/// @brief Aligned bounding-box helper kernels used by acceleration structures.
/// @details These methods intentionally keep vector mask support because they
/// operate on bounding-box batches in navigation and voxel managers, unlike the
/// scalar `UnplacedBox` navigation kernels.
struct ABBoxImplementation {

  /// @brief Intersect one ray with one scalar or lane-vector group of aligned boxes.
  /// @details The caller caches inverse direction and sign bits once per ray;
  /// vector callers use this to prune ABBox batches without touching the
  /// ordinary `UnplacedBox` navigation kernels.
  template <typename Real_v, typename basep>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static Real_v IntersectCachedKernel2(Vector3D<Real_v> const *corners,
                                                                                    Vector3D<basep> const &point,
                                                                                    Vector3D<basep> const &inverseray,
                                                                                    int signx, int signy, int signz,
                                                                                    basep t0, basep t1)
  {

    using Bool_v = vecCore::Mask<Real_v>;

    Real_v tmin  = (corners[signx].x() - Real_v(point.x())) * inverseray.x();
    Real_v tymax = (corners[1 - signy].y() - Real_v(point.y())) * inverseray.y();
    Bool_v done  = tmin > tymax;
    if (vecCore::MaskFull(done)) return InfinityLength<Real_v>();

    Real_v tmax  = (corners[1 - signx].x() - Real_v(point.x())) * inverseray.x();
    Real_v tymin = (corners[signy].y() - Real_v(point.y())) * inverseray.y();

    done |= (tymin > tmax);
    if (vecCore::MaskFull(done)) return InfinityLength<Real_v>();

    tmin = Max(tmin, tymin);
    tmax = Min(tmax, tymax);

    Real_v tzmin = (corners[signz].z() - point.z()) * inverseray.z();
    Real_v tzmax = (corners[1 - signz].z() - point.z()) * inverseray.z();

    done |= (Real_v(tmin) > Real_v(tzmax)) || (Real_v(tzmin) > Real_v(tmax));
    if (vecCore::MaskFull(done)) return InfinityLength<Real_v>();

    tmin = Max(tmin, tzmin);
    tmax = Min(tmax, tzmax);

    done |= !((tmin <= Real_v(t1 + kTolerance)) && (tmax > Real_v(t0 - kTolerance)));
    vecCore__MaskedAssignFunc(tmin, done, InfinityLength<Real_v>());
    return tmin;
  }

  // a contains kernel to be used with aligned bounding boxes
  // scalar and vector modes (aka backend) for boxes but only single points
  // should be useful to test one point against many bounding boxes
  // TODO: check if this can be unified with the normal generic box kernel
  template <typename Real_v, typename Bool_v = vecCore::Mask<Real_v>>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void ABBoxContainsKernel(Vector3D<Real_v> const &lowercorner,
                                                                               Vector3D<Real_v> const &uppercorner,
                                                                               Vector3D<Precision> const &point,
                                                                               Bool_v &inside)
  {

    inside = lowercorner.x() < Real_v(point.x());
    inside &= uppercorner.x() > Real_v(point.x());
    if (vecCore::MaskEmpty(inside)) return;

    inside &= lowercorner.y() < Real_v(point.y());
    inside &= uppercorner.y() > Real_v(point.y());
    if (vecCore::MaskEmpty(inside)) return;

    inside &= lowercorner.z() < Real_v(point.z());
    inside &= uppercorner.z() > Real_v(point.z());
  }

  // playing with a kernel that can do multi-box - single particle; multi-box -- multi-particle, single-box --
  // multi-particle
  template <typename T1, typename T2, typename Bool_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void ABBoxContainsKernelGeneric(Vector3D<T1> const &lowercorner,
                                                                                      Vector3D<T1> const &uppercorner,
                                                                                      Vector3D<T2> const &point,
                                                                                      Bool_v &inside)
  {
    inside = lowercorner.x() < T1(point.x());
    inside &= uppercorner.x() > T1(point.x());
    if (vecCore::MaskEmpty(inside)) return;

    inside &= lowercorner.y() < T1(point.y());
    inside &= uppercorner.y() > T1(point.y());
    if (vecCore::MaskEmpty(inside)) return;

    inside &= lowercorner.z() < T1(point.z());
    inside &= uppercorner.z() > T1(point.z());
  }

  // safety square for Bounding boxes
  // generic kernel treating one track and one or multiple boxes
  // in case a point is inside a box a squared value
  // is returned but given an overall negative sign
  template <typename Real_v, typename Real_s = typename vecCore::TypeTraits<Real_v>::ScalarType>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static Real_v ABBoxSafetySqr(Vector3D<Real_v> const &lowercorner,
                                                                            Vector3D<Real_v> const &uppercorner,
                                                                            Vector3D<Real_s> const &point)
  {

    using Vector3D_v = Vector3D<Real_v>;
    using Bool_v     = vecCore::Mask<Real_v>;

    const Vector3D_v kHalf(Real_v(static_cast<Real_s>(0.5)));
    const Vector3D_v origin((uppercorner + lowercorner) * kHalf);
    const Vector3D_v delta((uppercorner - lowercorner) * kHalf);
    // promote scalar point to vector point
    const Vector3D_v promotedpoint(Real_v(point.x()), Real_v(point.y()), Real_v(point.z()));

    // it would be nicer to have a standalone Abs function taking Vector3D as input
    const Vector3D_v safety = ((promotedpoint - origin).Abs()) - delta;
    const Bool_v outsidex   = safety.x() > Real_s(0.);
    const Bool_v outsidey   = safety.y() > Real_s(0.);
    const Bool_v outsidez   = safety.z() > Real_s(0.);

    Real_v runningsafetysqr(0.);                  // safety squared from outside
    Real_v runningmax(-InfinityLength<Real_v>()); // relevant for safety when we are inside

    // loop over dimensions manually unrolled
    // treat x dim
    {
      // this will be much simplified with operator notation
      Real_v tmp(0.);
      vecCore__MaskedAssignFunc(tmp, outsidex, safety.x() * safety.x());
      runningsafetysqr += tmp;
      runningmax = Max(runningmax, safety.x());
    }

    // treat y dim
    {
      Real_v tmp(0.);
      vecCore__MaskedAssignFunc(tmp, outsidey, safety.y() * safety.y());
      runningsafetysqr += tmp;
      runningmax = Max(runningmax, safety.y());
    }

    // treat z dim
    {
      Real_v tmp(0.);
      vecCore__MaskedAssignFunc(tmp, outsidez, safety.z() * safety.z());
      runningsafetysqr += tmp;
      runningmax = Max(runningmax, safety.z());
    }

    Bool_v inside = !(outsidex || outsidey || outsidez);
    if (!vecCore::MaskEmpty(inside)) vecCore__MaskedAssignFunc(runningsafetysqr, inside, -runningmax * runningmax);
    return runningsafetysqr;
  }

  // safety square for Bounding boxes, returning the squared range for any point in the box
  // generic kernel treating one track and one or multiple boxes
  // in case a point is inside a box a squared value
  // is returned but given an overall negative sign
  template <typename Real_v, typename Real_s = typename vecCore::TypeTraits<Real_v>::ScalarType>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static Real_v ABBoxSafetyRangeSqr(Vector3D<Real_v> const &lowercorner,
                                                                                 Vector3D<Real_v> const &uppercorner,
                                                                                 Vector3D<Real_s> const &point,
                                                                                 Real_v &safetymaxsqr)
  {

    using Vector3D_v = Vector3D<Real_v>;
    using Bool_v     = vecCore::Mask<Real_v>;

    const Vector3D_v kHalf(Real_v(static_cast<Real_s>(0.5)));
    const Vector3D_v origin((uppercorner + lowercorner) * kHalf);
    const Vector3D_v delta((uppercorner - lowercorner) * kHalf);
    // promote scalar point to vector point
    const Vector3D_v promotedpoint(Real_v(point.x()), Real_v(point.y()), Real_v(point.z()));

    // it would be nicer to have a standalone Abs function taking Vector3D as input
    const Vector3D_v safety  = ((promotedpoint - origin).Abs()) - delta;
    const Vector3D_v safetyp = ((promotedpoint - origin).Abs()) + delta;
    const Bool_v outsidex    = safety.x() > Real_s(0.);
    const Bool_v outsidey    = safety.y() > Real_s(0.);
    const Bool_v outsidez    = safety.z() > Real_s(0.);

    Real_v runningsafetysqr(0.);                  // safety squared from outside
    safetymaxsqr = safetyp.Mag2();                // safetymax squared from outside
    Real_v runningmax(-InfinityLength<Real_v>()); // relevant for safety when we are inside

    // loop over dimensions manually unrolled
    // treat x dim
    {
      // this will be much simplified with operator notation
      Real_v tmp(0.);
      vecCore__MaskedAssignFunc(tmp, outsidex, safety.x() * safety.x());
      runningsafetysqr += tmp;
      runningmax = Max(runningmax, safety.x());
      //      vecCore__MaskedAssignFunc(tmp, outsidex, safetyp.x() * safetyp.x());
      //      safetymaxsqr += tmp;
    }

    // treat y dim
    {
      Real_v tmp(0.);
      vecCore__MaskedAssignFunc(tmp, outsidey, safety.y() * safety.y());
      runningsafetysqr += tmp;
      runningmax = Max(runningmax, safety.y());
      //      vecCore__MaskedAssignFunc(tmp, outsidey, safetyp.y() * safetyp.y());
      //      safetymaxsqr += tmp;
    }

    // treat z dim
    {
      Real_v tmp(0.);
      vecCore__MaskedAssignFunc(tmp, outsidez, safety.z() * safety.z());
      runningsafetysqr += tmp;
      runningmax = Max(runningmax, safety.z());
      //      vecCore__MaskedAssignFunc(tmp, outsidez, safetyp.z() * safetyp.z());
      //      safetymaxsqr += tmp;
    }

    Bool_v inside = !(outsidex || outsidey || outsidez);
    if (!vecCore::MaskEmpty(inside)) vecCore__MaskedAssignFunc(runningsafetysqr, inside, -runningmax * runningmax);
    return runningsafetysqr;
  }

}; // end aligned bounding box struct

} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif // VECGEOM_MANAGEMENT_ABBOXIMPLEMENTATION_H_
