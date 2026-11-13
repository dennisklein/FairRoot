//===-- kernel/BoxImplementation.h ----------------------------------*- C++ -*-===//
//===--------------------------------------------------------------------------===//
/// @file BoxImplementation.h
/// @brief Axis-aligned box kernel helpers and navigation entry points.
/// @author Johannes de Fine Licht (johannes.definelicht@cern.ch), Sandro Wenzel (sandro.wenzel@cern.ch)

/// History notes:
/// 2013 - 2014: original development (abstracted kernels); Johannes and Sandro
/// Oct 2015: revision + moving to new backend structure (Sandro Wenzel)

#ifndef VECGEOM_VOLUMES_KERNEL_BOXIMPLEMENTATION_H_
#define VECGEOM_VOLUMES_KERNEL_BOXIMPLEMENTATION_H_

#include "VecGeom/base/Vector3D.h"
#include "VecGeom/volumes/BoxStruct.h"
#include "VecGeom/volumes/SurfaceHitView.h"
#include "VecGeom/volumes/kernel/GenericKernels.h"
#include <VecCore/VecCore>

#include <cstdio>

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE(struct BoxImplementation;);
VECGEOM_DEVICE_DECLARE_CONV(struct, BoxImplementation);

inline namespace VECGEOM_IMPL_NAMESPACE {

class PlacedBox;
template <typename T>
struct BoxStruct;
class UnplacedBox;

/// @brief Kernel implementation for axis-aligned boxes.
/// @details The stored dimensions are half-lengths along x, y, and z. The main
/// `UnplacedBox` navigation entry points are scalar and use local slab
/// predicates for classification, distance, safety, and normals.
struct BoxImplementation {

  using PlacedShape_t    = PlacedBox;
  using UnplacedStruct_t = BoxStruct<Precision>;
  using UnplacedVolume_t = UnplacedBox;

  /// @brief Return the box half-lengths promoted to the requested real type.
  /// @param box Box runtime data.
  /// @return Vector of x/y/z half-lengths.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static Vector3D<Real_v> HalfSize(const UnplacedStruct_t &box)
  {
    return Vector3D<Real_v>(box.fDimensions[0], box.fDimensions[1], box.fDimensions[2]);
  }

  /// @brief Test whether a point is contained by the box.
  /// @details The containment convention accepts surface points within the box
  /// tolerance, matching the other shape kernels.
  /// @param box Box runtime data.
  /// @param point Query point.
  /// @param[out] inside True when @p point is inside or on the tolerated surface.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void Contains(UnplacedStruct_t const &box,
                                                                    Vector3D<Real_v> const &point, bool &inside)
  {
    // in analogy to other shapes, surface points are considered inside
    inside = (point.Abs() - HalfSize<Real_v>(box)).Max() < Real_v(kTolerance);
  }

  /// @brief Classify a point as inside, outside, or on the box surface.
  /// @details Classification is based on the maximum signed distance to the
  /// three slabs. Points within `kHalfTolerance` of any limiting slab are
  /// reported as `kSurface`.
  /// @param box Box runtime data.
  /// @param point Query point.
  /// @param[out] inside VecGeom inside code.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void Inside(UnplacedStruct_t const &box,
                                                                  Vector3D<Real_v> const &point, Inside_t &inside)
  {
    Real_v dist = (point.Abs() - HalfSize<Real_v>(box)).Max();

    inside = (Abs(dist) < Real_v(kHalfTolerance)) ? Inside_t(kSurface)
                                                  : ((dist < Real_v(0.0)) ? Inside_t(kInside) : Inside_t(kOutside));
  }

  /// @brief Shared scalar slab classification helper.
  /// @details The maximum signed slab distance determines outside rejection. If
  /// `ForInside` is true, the helper also reports strict interior ownership
  /// separated from all slabs by at least `kHalfTolerance`.
  /// @param halfsize Box half-lengths.
  /// @param point Query point.
  /// @param[out] completelyinside True when the point is strictly inside.
  /// @param[out] completelyoutside True when the point is outside.
  template <typename Real_v, bool ForInside>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void GenericKernelForContainsAndInside(
      Vector3D<Real_v> const &halfsize, Vector3D<Real_v> const &point, bool &completelyinside, bool &completelyoutside)
  {
    Real_v dist = (point.Abs() - halfsize).Max();

    if (ForInside) completelyinside = dist < Real_v(-kHalfTolerance);

    completelyoutside = dist > Real_v(kHalfTolerance);
  }

  /// @brief Compute distance from an exterior point to the box.
  /// @details Uses the standard slab intersection interval. The entry distance
  /// is rejected when the interval is empty, when the candidate exit is already
  /// behind or on the start point, or when the start is already on an exit slab
  /// within tolerance.
  /// @param box Box runtime data.
  /// @param point Query point, expected outside.
  /// @param direction Normalized propagation direction.
  /// @param stepMax Unused by this implementation.
  /// @param[out] distance Entry distance, or infinity when no valid entry exists.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void DistanceToIn(UnplacedStruct_t const &box,
                                                                        Vector3D<Real_v> const &point,
                                                                        Vector3D<Real_v> const &direction,
                                                                        Real_v const & /* stepMax */, Real_v &distance)
  {
    const Vector3D<Real_v> invDir(Real_v(1.0) / NonZero(direction[0]), Real_v(1.0) / NonZero(direction[1]),
                                  Real_v(1.0) / NonZero(direction[2]));

    const Vector3D<Real_v> signDir(Sign(direction[0]), Sign(direction[1]), Sign(direction[2]));

    const Vector3D<Real_v> tempIn  = -signDir * box.fDimensions - point;
    const Vector3D<Real_v> tempOut = signDir * box.fDimensions - point;

    // add a check for point on exit surface
    const Real_v absOrthogOut = Abs((signDir * tempOut).Min());

    const Real_v distOut = (tempOut * invDir).Min();

    // distIn calculation
    distance = (tempIn * invDir).Max();

    if (distance >= distOut || distOut <= Real_v(kHalfTolerance) || absOrthogOut <= Real_v(kHalfTolerance)) {
      distance = InfinityLength<Real_v>();
    }
  }

  /// @brief Compute distance from an interior point to leave the box and optionally its face normal.
  /// @details The nearest forward slab crossing is selected after skipping axes
  /// whose direction is effectively parallel for the current half-length scale.
  /// Points already outside the tolerated box return `-1`. A requested scalar
  /// normal reuses the winning slab distances without reclassifying the hit point.
  /// @param box Box runtime data.
  /// @param point Query point, expected inside or on the surface.
  /// @param direction Normalized propagation direction.
  /// @param stepMax Unused by this implementation.
  /// @param[out] distance Exit distance, or `-1` for outside input.
  /// @param hit_info Optional scalar exit-normal destination and supporting-plane flag.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void DistanceToOut(UnplacedStruct_t const &box,
                                                                         Vector3D<Real_v> const &point,
                                                                         Vector3D<Real_v> const &direction,
                                                                         Real_v const & /* stepMax */, Real_v &distance,
                                                                         SurfaceHitView<Real_v> *hit_info = nullptr)
  {
    if (hit_info) hit_info->Clear();
    distance = Real_v(-1.0);

    // Quick check: if the point is outside the box plus tolerance, return early
    const Real_v pointDistance = (point.Abs() - HalfSize<Real_v>(box)).Max();
    if (pointDistance > Real_v(kTolerance)) return;

    // Compute one complete slab expression at a time, keeping the reciprocal form.
    const Vector3D<Real_v> slabDistance(
        (Sign(direction[0]) * box.fDimensions[0] - point[0]) * (Real_v(1.) / NonZero(direction[0])),
        (Sign(direction[1]) * box.fDimensions[1] - point[1]) * (Real_v(1.) / NonZero(direction[1])),
        (Sign(direction[2]) * box.fDimensions[2] - point[2]) * (Real_v(1.) / NonZero(direction[2])));

    // Skip near-parallel directions, whose exit distances are numerically unreliable.
    auto skip = direction.Abs() * box.fDimensions < Vector3D<Real_v>(kTolerance);

    // Compute final distance using min of unskipped axes
    distance = slabDistance.MinSkip(skip);

    // A tolerated surface start can produce a small negative slab distance.
    // Follow Inside() for hit validity instead of rejecting by distance sign.
    if (hit_info && pointDistance < Real_v(kHalfTolerance)) {
      if (hit_info->WantsNormal()) {
        Vector3D<Real_v> normal(Real_v(0.));
        if (!skip[2] && slabDistance[2] == distance) {
          normal[2] = Sign(direction[2]);
        } else if (!skip[0] && slabDistance[0] == distance) {
          normal[0] = Sign(direction[0]);
        } else {
          normal[1] = Sign(direction[1]);
        }
        hit_info->SetNormal(normal);
      }
      hit_info->SetLeftBehind(true);
    }
  }

  /// @brief Compute the safety from an exterior point to the box.
  /// @details The value is the maximum signed slab excess and is non-positive
  /// for points already inside the box.
  /// @param box Box runtime data.
  /// @param point Query point.
  /// @param[out] safety Conservative distance to the nearest entry boundary.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SafetyToIn(UnplacedStruct_t const &box,
                                                                      Vector3D<Real_v> const &point, Real_v &safety)
  {
    safety = (point.Abs() - HalfSize<Real_v>(box)).Max();
  }

  /// @brief Compute the safety from an interior point to leave the box.
  /// @details The value is the minimum remaining distance to any of the six
  /// limiting slabs and becomes negative for outside points.
  /// @param box Box runtime data.
  /// @param point Query point.
  /// @param[out] safety Conservative distance to the nearest exit boundary.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SafetyToOut(UnplacedStruct_t const &box,
                                                                       Vector3D<Real_v> const &point, Real_v &safety)
  {
    safety = (HalfSize<Real_v>(box) - point.Abs()).Min();
  }

  /// @brief Compute a box surface normal.
  /// @details A point within `kHalfTolerance` of one face returns that face
  /// normal. Edge and corner points average all equally closest face normals and
  /// normalize the result.
  /// @param box Box runtime data.
  /// @param point Query point.
  /// @param[out] valid True when @p point is close enough to a box surface.
  /// @return Outward normal. For invalid input this is a fallback vector.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static Vector3D<Real_v> NormalKernel(UnplacedStruct_t const &box,
                                                                                    Vector3D<Real_v> const &point,
                                                                                    bool &valid)
  {
    // Computes the normal on a surface and returns it as a unit vector
    //   In case a point is further than kHalfTolerance from a surface, set valid=false
    //   Must return a valid vector. (even if the point is not on the surface.)
    //
    //   On an edge or corner, provide an average normal of all facets within tolerance

    const Vector3D<Real_v> safety((point.Abs() - HalfSize<Real_v>(box)).Abs());
    const Real_v safmin = safety.Min();
    valid               = safmin < kHalfTolerance;

    Vector3D<Real_v> normal(0.);
    if (safety[0] - safmin < kHalfTolerance) normal[0] = Sign(point[0]);
    if (safety[1] - safmin < kHalfTolerance) normal[1] = Sign(point[1]);
    if (safety[2] - safmin < kHalfTolerance) normal[2] = Sign(point[2]);
    if (normal.Mag2() > 1.0) normal.Normalize();

    return normal;
  }

  /// @brief Slab-intersection helper for one scalar ray and one box.
  /// @details This legacy helper is used by assembly and bounding-box code. It
  /// tests whether the ray intersects the box described by the lower/upper
  /// corners passed in @p corners.
  /// @param corners Two box corners ordered as lower and upper.
  /// @param point Ray start point.
  /// @param ray Ray direction.
  /// @return True if the ray intersects the box interval.
  // an algorithm to test for intersection ( could be faster than DistanceToIn )
  // actually this also calculated the distance at the same time ( in tmin )
  // template <class Backend>
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  static bool Intersect(Vector3D<Precision> const *corners, Vector3D<Precision> const &point,
                        Vector3D<Precision> const &ray, Precision /* t0 */, Precision /* t1 */)
  {
    // intersection algorithm 1 ( Amy Williams )
    Precision tmin, tmax, tymin, tymax, tzmin, tzmax;

    // IF THERE IS A STEPMAX; COULD ALSO CHECK SAFETIES
    Precision inverserayx = 1. / ray[0];
    Precision inverserayy = 1. / ray[1];

    // TODO: we should promote this to handle multiple boxes
    int sign[3];
    sign[0] = inverserayx < 0;
    sign[1] = inverserayy < 0;

    tmin  = (corners[sign[0]].x() - point.x()) * inverserayx;
    tmax  = (corners[1 - sign[0]].x() - point.x()) * inverserayx;
    tymin = (corners[sign[1]].y() - point.y()) * inverserayy;
    tymax = (corners[1 - sign[1]].y() - point.y()) * inverserayy;

    if ((tmin > tymax) || (tymin > tmax)) return false;

    Precision inverserayz = 1. / ray.z();
    sign[2]               = inverserayz < 0;

    if (tymin > tmin) tmin = tymin;
    if (tymax < tmax) tmax = tymax;

    tzmin = (corners[sign[2]].z() - point.z()) * inverserayz;
    tzmax = (corners[1 - sign[2]].z() - point.z()) * inverserayz;

    if ((tmin > tzmax) || (tzmin > tmax)) return false;
    // if ((tzmin > tmin)) tmin = tzmin;
    // if (tzmax < tmax) tmax   = tzmax;
    // return ((tmin < t1) && (tmax > t0));
    // std::cerr << "tmin " << tmin << " tmax " << tmax << "\n";
    return true;
  }

  // an algorithm to test for intersection ( could be faster than DistanceToIn )
  // actually this also calculated the distance at the same time ( in tmin )
  template <int signx, int signy, int signz>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE
      //__attribute__((noinline))
      static Precision IntersectCached(Vector3D<Precision> const *corners, Vector3D<Precision> const &point,
                                       Vector3D<Precision> const &inverseray, Precision t0, Precision t1)
  {
    // intersection algorithm 1 ( Amy Williams )

    // NOTE THE FASTEST VERSION IS STILL THE ORIGINAL IMPLEMENTATION

    Precision tmin, tmax, tymin, tymax, tzmin, tzmax;

    // TODO: we should promote this to handle multiple boxes
    // observation: we always compute sign and 1-sign; so we could do the assignment
    // to tmin and tmax in a masked assignment thereafter
    tmin  = (corners[signx].x() - point.x()) * inverseray.x();
    tmax  = (corners[1 - signx].x() - point.x()) * inverseray.x();
    tymin = (corners[signy].y() - point.y()) * inverseray.y();
    tymax = (corners[1 - signy].y() - point.y()) * inverseray.y();
    if ((tmin > tymax) || (tymin > tmax)) return InfinityLength<Precision>();

    if (tymin > tmin) tmin = tymin;
    if (tymax < tmax) tmax = tymax;

    tzmin = (corners[signz].z() - point.z()) * inverseray.z();
    tzmax = (corners[1 - signz].z() - point.z()) * inverseray.z();

    if ((tmin > tzmax) || (tzmin > tmax)) return InfinityLength<Precision>(); // false
    if ((tzmin > tmin)) tmin = tzmin;
    if (tzmax < tmax) tmax = tzmax;

    if (!((tmin < t1) && (tmax > t0))) return InfinityLength<Precision>();
    return tmin;
  }

}; // End struct BoxImplementation
} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif // VECGEOM_VOLUMES_KERNEL_BOXIMPLEMENTATION_H_
