// This file is part of VecGeom and is distributed under the
// conditions in the file LICENSE.txt in the top directory.
// For the full list of authors see CONTRIBUTORS.txt and `git log`.

/// @file volumes/kernel/EllipsoidImplementation.h
/// @brief Navigation kernels for the ellipsoid solid.
/// @author Evgueni Tcherniaev

#ifndef VECGEOM_VOLUMES_KERNEL_ELLIPSOIDIMPLEMENTATION_H_
#define VECGEOM_VOLUMES_KERNEL_ELLIPSOIDIMPLEMENTATION_H_

#include "VecGeom/base/Vector3D.h"
#include "VecGeom/volumes/EllipsoidStruct.h"
#include "VecGeom/volumes/SurfaceHitView.h"
#include "VecGeom/volumes/kernel/GenericKernels.h"
#include <VecCore/VecCore>

#include <cstdio>
#include <iomanip>

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE(struct EllipsoidImplementation;);
VECGEOM_DEVICE_DECLARE_CONV(struct, EllipsoidImplementation);

inline namespace VECGEOM_IMPL_NAMESPACE {

class PlacedEllipsoid;
template <typename T>
struct EllipsoidStruct;
class UnplacedEllipsoid;

/// @brief Implements Ellipsoid classification, distance, safety, and normal kernels.
///
/// @details The implementation maps the ellipsoid to a sphere using cached scale
/// factors from `EllipsoidStruct`. Z cuts are handled in the scaled frame for
/// classification and distance root filtering, while safety values are converted
/// back to conservative local-space bounds where needed.
struct EllipsoidImplementation {

  using PlacedShape_t    = PlacedEllipsoid;
  using UnplacedStruct_t = EllipsoidStruct<Precision>;
  using UnplacedVolume_t = UnplacedEllipsoid;

  /// @brief Store the supporting normal for a selected ellipsoid exit.
  /// @tparam Real_v Scalar floating-point type used by the distance kernel.
  /// @param ellipsoid Cached unplaced ellipsoid data.
  /// @param px Scaled start x coordinate.
  /// @param py Scaled start y coordinate.
  /// @param pz Scaled start z coordinate.
  /// @param vx Scaled direction x component.
  /// @param vy Scaled direction y component.
  /// @param vz Scaled direction z component.
  /// @param distance Returned exit distance.
  /// @param hitZ True when a z-cut plane supplied the exit.
  /// @param[out] hitInfo Non-null requested normal destination and supporting-plane result.
  /// @details A clipped ellipsoid is convex, so both cut-plane normals and
  /// curved-surface gradient normals define supporting exit planes.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SetExitHit(UnplacedStruct_t const &ellipsoid, Real_v px,
                                                                      Real_v py, Real_v pz, Real_v vx, Real_v vy,
                                                                      Real_v vz, Real_v distance, bool hitZ,
                                                                      SurfaceHitView<Real_v> *hitInfo)
  {
    if (hitInfo->WantsNormal()) {
      Vector3D<Real_v> normal;
      if (hitZ) {
        normal.Set(Real_v(0.), Real_v(0.), Sign(pz + distance * vz - ellipsoid.fScZMidCut));
      } else {
        normal.Set((px + distance * vx) * ellipsoid.fSx, (py + distance * vy) * ellipsoid.fSy,
                   (pz + distance * vz) * ellipsoid.fSz);
        normal.Normalize();
      }
      hitInfo->SetNormal(normal);
    }
    hitInfo->SetLeftBehind(true);
  }

  /// @brief Tests whether a point is inside or on the ellipsoid.
  /// @param ellipsoid Cached unplaced ellipsoid data.
  /// @param point Point in the local solid frame.
  /// @param inside Output flag set to true unless the point is outside.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void Contains(UnplacedStruct_t const &ellipsoid,
                                                                    Vector3D<Real_v> const &point, bool &inside)
  {
    bool unused(false), outside(false);
    GenericKernelForContainsAndInside<Real_v, false>(ellipsoid, point, unused, outside);
    inside = !outside;
  }

  /// @brief Classifies a point against the ellipsoid and its z cuts.
  /// @param ellipsoid Cached unplaced ellipsoid data.
  /// @param point Point in the local solid frame.
  /// @param inside Output classification as `EInside`.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void Inside(UnplacedStruct_t const &ellipsoid,
                                                                  Vector3D<Real_v> const &point, Inside_t &inside)
  {
    bool completelyinside(false), completelyoutside(false);
    GenericKernelForContainsAndInside<Real_v, true>(ellipsoid, point, completelyinside, completelyoutside);
    inside = EInside::kSurface;
    if (completelyoutside) inside = Inside_t(EInside::kOutside);
    if (completelyinside) inside = Inside_t(EInside::kInside);
  }

  /// @brief Computes the shared inside/outside predicates.
  /// @param ellipsoid Cached unplaced ellipsoid data.
  /// @param point Point in the local solid frame.
  /// @param completelyinside Output true when the point is strictly inside all bounds.
  /// @param completelyoutside Output true when the point is outside at least one bound.
  /// @tparam ForInside Whether the strict-inside predicate is needed.
  ///
  /// @details The point is scaled to the auxiliary sphere. The radial and z-cut
  /// distances are combined into a single tolerance-band predicate.
  template <typename Real_v, bool ForInside>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void GenericKernelForContainsAndInside(
      UnplacedStruct_t const &ellipsoid, Vector3D<Real_v> const &point, bool &completelyinside, bool &completelyoutside)
  {
    Real_v x      = point.x() * ellipsoid.fSx;
    Real_v y      = point.y() * ellipsoid.fSy;
    Real_v z      = point.z() * ellipsoid.fSz;
    Real_v distZ  = Abs(z - ellipsoid.fScZMidCut) - ellipsoid.fScZDimCut;
    Real_v distR  = ellipsoid.fQ1 * (x * x + y * y + z * z) - ellipsoid.fQ2;
    Real_v safety = Max(distZ, distR);

    completelyoutside = safety > kHalfTolerance;
    if (ForInside) completelyinside = safety <= -kHalfTolerance;
    return;
  }

  /// @brief Computes the distance from an outside point to the ellipsoid.
  /// @param ellipsoid Cached unplaced ellipsoid data.
  /// @param point Point in the local solid frame.
  /// @param direction Unit direction in the local solid frame.
  /// @param distance Output distance, or `kInfLength` if the ray misses.
  ///
  /// @details Far-away points are translated closer along the ray before scaling
  /// to the auxiliary sphere. Candidate hits are the overlap of the sphere
  /// interval and the z-cut interval.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void DistanceToIn(UnplacedStruct_t const &ellipsoid,
                                                                        Vector3D<Real_v> const &point,
                                                                        Vector3D<Real_v> const &direction,
                                                                        Real_v const & /*stepMax*/, Real_v &distance)
  {
    distance = kInfLength;
    Real_v offset(0.);
    Vector3D<Real_v> pcur(point);

    // Move point closer, if required
    Real_v Rfar2(1024. * ellipsoid.fRsph * ellipsoid.fRsph); // 1024 = 32 * 32
    if ((pcur.Mag2() > Rfar2) && (direction.Dot(point) < Real_v(0.))) {
      offset = pcur.Mag() - Real_v(2.) * ellipsoid.fRsph;
      pcur += offset * direction;
    }

    // Scale ellipsoid to sphere
    Real_v px = pcur.x() * ellipsoid.fSx;
    Real_v py = pcur.y() * ellipsoid.fSy;
    Real_v pz = pcur.z() * ellipsoid.fSz;
    Real_v vx = direction.x() * ellipsoid.fSx;
    Real_v vy = direction.y() * ellipsoid.fSy;
    Real_v vz = direction.z() * ellipsoid.fSz;

    // Check if point is leaving the solid
    Real_v pzcut = pz - ellipsoid.fScZMidCut;
    Real_v dzcut = Real_v(ellipsoid.fScZDimCut);
    Real_v distZ = Abs(pzcut) - dzcut;

    Real_v rr    = px * px + py * py + pz * pz;
    Real_v vv    = vx * vx + vy * vy + vz * vz;
    Real_v pv    = px * vx + py * vy + pz * vz;
    Real_v distR = ellipsoid.fQ1 * rr - ellipsoid.fQ2;
    bool leaving =
        (distZ >= -kHalfTolerance && pzcut * vz >= Real_v(0.)) || (distR >= -kHalfTolerance && pv >= Real_v(0.));
    if (leaving) return;

    // Find intersection with Z planes
    Real_v invz  = Real_v(-1.) / NonZero(vz);
    Real_v dz    = CopySign(dzcut, invz);
    Real_v tzmin = (pzcut - dz) * invz;
    Real_v tzmax = (pzcut + dz) * invz;

    // Find intersection with sphere
    Real_v A   = vv;
    Real_v B   = pv;
    Real_v C   = (rr - ellipsoid.fR * ellipsoid.fR);
    Real_v D   = B * B - A * C;
    Real_v EPS = Real_v(2.) * rr * vv * kEpsilon;
    if (D <= EPS) {
      if (C <= Real_v(0.) && distR > -kTolerance && distR < kTolerance && B < -Real_v(0.5 * kTolerance) * A)
        distance = offset;
      return;
    }
    Real_v sqrtD = Sqrt(D);
    Real_v trmin = (-B - sqrtD) / A;
    Real_v trmax = (-B + sqrtD) / A;

    // Set preliminary distances to in/out
    Real_v tmin = Max(tzmin, trmin);
    Real_v tmax = Min(tzmax, trmax);

    // Check if no intersection
    bool done = (tmax - tmin) <= kTolerance;

    // Set distance
    if (!done) distance = tmin + offset;
  }

  /// @brief Computes the distance from an inside point to the next boundary.
  /// @tparam Real_v Scalar floating-point type used by the distance kernel.
  /// @param ellipsoid Cached unplaced ellipsoid data.
  /// @param point Point in the local solid frame.
  /// @param direction Unit direction in the local solid frame.
  /// @param stepMax Unused by this implementation.
  /// @param[out] distance Exit distance, zero for unresolved tangent/no-hit cases, or -1 for wrong-side input.
  /// @param[out] hitInfo Optional exit-normal destination and supporting-plane result.
  ///
  /// @details The first exit is the smaller positive hit between the scaled
  /// sphere and the active z-cut plane. A requested normal reuses that winning
  /// candidate and computes the curved gradient only when needed.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void DistanceToOut(UnplacedStruct_t const &ellipsoid,
                                                                         Vector3D<Real_v> const &point,
                                                                         Vector3D<Real_v> const &direction,
                                                                         Real_v const & /* stepMax */, Real_v &distance,
                                                                         SurfaceHitView<Real_v> *hitInfo = nullptr)
  {
    if (hitInfo) hitInfo->Clear();

    // Scale ellipsoid to sphere
    Real_v px = point.x() * ellipsoid.fSx;
    Real_v py = point.y() * ellipsoid.fSy;
    Real_v pz = point.z() * ellipsoid.fSz;
    Real_v vx = direction.x() * ellipsoid.fSx;
    Real_v vy = direction.y() * ellipsoid.fSy;
    Real_v vz = direction.z() * ellipsoid.fSz;

    // Check if point is outside ("wrong side")
    Real_v pzcut = pz - ellipsoid.fScZMidCut;
    Real_v dzcut = Real_v(ellipsoid.fScZDimCut);
    Real_v distZ = Abs(pzcut) - dzcut;

    Real_v rr    = px * px + py * py + pz * pz;
    Real_v vv    = vx * vx + vy * vy + vz * vz;
    Real_v pv    = px * vx + py * vy + pz * vz;
    Real_v distR = ellipsoid.fQ1 * rr - ellipsoid.fQ2;
    bool outside = Max(distR, distZ) > kHalfTolerance;

    distance = Real_v(0.);
    if (outside) {
      distance = Real_v(-1.);
      return;
    }

    // Find intersection with Z planes
    Real_v tzmax = kMaximum;
    if (vz != Real_v(0.)) tzmax = (CopySign(Real_v(dzcut), vz) - pzcut) / vz;

    // Find intersection with sphere
    Real_v B   = pv / vv;
    Real_v C   = (rr - ellipsoid.fR * ellipsoid.fR) / vv;
    Real_v D   = B * B - C;
    Real_v EPS = Real_v(2.) * rr * vv * kEpsilon;
    if (D <= EPS) {
      const bool hasCurvedExit =
          C <= Real_v(0.) && distR > -kTolerance && distR < kTolerance && B < -Real_v(0.5 * kTolerance);
      if (hasCurvedExit) {
        const Real_v curvedDistance = -Real_v(2.) * B;
        distance                    = Min(tzmax, curvedDistance);
        if (hitInfo) SetExitHit(ellipsoid, px, py, pz, vx, vy, vz, distance, tzmax <= curvedDistance, hitInfo);
      } else if (hitInfo) {
        const bool onCurved = distR > -kTolerance && distR < kTolerance;
        const bool onZ      = distZ > -kTolerance && distZ < kTolerance;
        if (onCurved || onZ) SetExitHit(ellipsoid, px, py, pz, vx, vy, vz, distance, !onCurved && onZ, hitInfo);
      }
      return;
    }
    Real_v sqrtD = Sqrt(D);
    Real_v trmax = -B + sqrtD;

    // Set distance
    distance = Min(tzmax, trmax);
    if (hitInfo) SetExitHit(ellipsoid, px, py, pz, vx, vy, vz, distance, tzmax <= trmax, hitInfo);
  }

  /// @brief Computes a conservative safety from an outside point to the ellipsoid.
  /// @param ellipsoid Cached unplaced ellipsoid data.
  /// @param point Point in the local solid frame.
  /// @param safety Output lower bound to the solid boundary.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SafetyToIn(UnplacedStruct_t const &ellipsoid,
                                                                      Vector3D<Real_v> const &point, Real_v &safety)
  {
    Real_v x = point.x() * ellipsoid.fSx;
    Real_v y = point.y() * ellipsoid.fSy;
    Real_v z = point.z() * ellipsoid.fSz;
    Real_v r = Sqrt(x * x + y * y + z * z);
    // Set safety to zero if point is on surface
    Real_v safeZ = Abs(z - ellipsoid.fScZMidCut) - ellipsoid.fScZDimCut;
    Real_v safeR = r - ellipsoid.fR;
    safety       = Max(safeZ, safeR);
    if (Abs(safety) <= kHalfTolerance) safety = Real_v(0.);
    // Adjust safety using bounding box
    Real_v distZ  = Max(point.z() - ellipsoid.fZTopCut, ellipsoid.fZBottomCut - point.z());
    Real_v distXY = Max(Abs(point.x()) - ellipsoid.fXmax, Abs(point.y()) - ellipsoid.fYmax);
    if (safety > Real_v(0.)) safety = Max(safety, distZ, distXY);
  }

  /// @brief Computes a conservative safety from an inside point to the boundary.
  /// @param ellipsoid Cached unplaced ellipsoid data.
  /// @param point Point in the local solid frame.
  /// @param safety Output lower bound to the exit boundary.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SafetyToOut(UnplacedStruct_t const &ellipsoid,
                                                                       Vector3D<Real_v> const &point, Real_v &safety)
  {
    Real_v x = point.x() * ellipsoid.fSx;
    Real_v y = point.y() * ellipsoid.fSy;
    Real_v z = point.z() * ellipsoid.fSz;
    // Set safety to zero if point is on surface
    Real_v safeR = ellipsoid.fR - Sqrt(x * x + y * y + z * z);
    Real_v safeZ = ellipsoid.fScZDimCut - Abs(z - ellipsoid.fScZMidCut);
    safety       = Min(safeZ, safeR);
    if (Abs(safety) <= kHalfTolerance) {
      safety = Real_v(0.);
      return;
    }
    // Only positive safety needs the unscaled cap correction.
    if (safety > Real_v(0.)) {
      Real_v distZ = Min(ellipsoid.fZTopCut - point.z(), point.z() - ellipsoid.fZBottomCut);
      safety       = Min(safeR, distZ);
    }
  }

  /// @brief Computes the outward normal on an ellipsoid boundary point.
  /// @param ellipsoid Cached unplaced ellipsoid data.
  /// @param point Point in the local solid frame.
  /// @param valid Output false when the point is not within tolerance of a boundary.
  /// @return Unit normal for boundary points, or a nearest-boundary fallback normal.
  ///
  /// @details Side normals use the ellipsoid gradient. Points on a z-cut edge
  /// receive the normalized average of side and cut normals.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static Vector3D<Real_v> NormalKernel(UnplacedStruct_t const &ellipsoid,
                                                                                    Vector3D<Real_v> const &point,
                                                                                    bool &valid)
  {
    Vector3D<Real_v> normal(0.);
    valid = true;

    Real_v px   = point.x();
    Real_v py   = point.y();
    Real_v pz   = point.z();
    Real_v x    = px * ellipsoid.fSx;
    Real_v y    = py * ellipsoid.fSy;
    Real_v z    = pz * ellipsoid.fSz;
    Real_v mag2 = x * x + y * y + z * z;

    // The scaled gradient differs by the positive factor fR squared, removed by Unit().
    // Check lateral surface
    Real_v distR      = ellipsoid.fQ1 * mag2 - ellipsoid.fQ2;
    const bool onSide = Abs(distR) <= kHalfTolerance;
    if (onSide) {
      normal = Vector3D<Real_v>(x * ellipsoid.fSx, y * ellipsoid.fSy, z * ellipsoid.fSz).Unit();
    }

    // Check z cuts
    Real_v distZ = Abs(z - ellipsoid.fScZMidCut) - ellipsoid.fScZDimCut;
    if (Abs(distZ) <= kHalfTolerance) {
      normal[2] += Sign(z - ellipsoid.fScZMidCut);
    } else if (onSide) {
      // The side gradient is already normalized; only a cut edge may need averaging.
      return normal;
    }

    // Average normal, if required
    Real_v normalMag2 = normal.Mag2();
    if (normalMag2 > 1.) {
      normal     = normal.Unit();
      normalMag2 = Real_v(1.);
    }
    if (normalMag2 > Real_v(0.)) return normal;

    // Point is not on the surface - normally, this should never be
    // Return normal to the nearest surface
    valid     = false;
    normal[2] = Sign(z - ellipsoid.fScZMidCut);
    distR     = Sqrt(mag2) - ellipsoid.fR;
    if (distR > distZ && mag2 > Real_v(0.)) {
      normal = Vector3D<Real_v>(x * ellipsoid.fSx, y * ellipsoid.fSy, z * ellipsoid.fSz).Unit();
    }
    return normal;
  }
};
} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif // VECGEOM_VOLUMES_KERNEL_ELLIPSOIDIMPLEMENTATION_H_
