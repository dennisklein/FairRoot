// This file is part of VecGeom and is distributed under the
// conditions in the file LICENSE.txt in the top directory.
// For the full list of authors see CONTRIBUTORS.txt and `git log`.

/// @file OrbImplementation.h
/// @brief Navigation kernels for the Orb solid.
/// \author Raman Sehgal

/// History notes:
/// 2014 - 2015: original development (abstracted kernels); Raman Sehgal
/// July 2016: revision + moving to new backend structure (Raman Sehgal)

#ifndef VECGEOM_VOLUMES_KERNEL_ORBIMPLEMENTATION_H_
#define VECGEOM_VOLUMES_KERNEL_ORBIMPLEMENTATION_H_

#include "VecGeom/base/Vector3D.h"
#include "VecGeom/volumes/OrbStruct.h"
#include "VecGeom/volumes/SurfaceHitView.h"
#include "VecGeom/volumes/kernel/GenericKernels.h"
#include <VecCore/VecCore>

#include <cstdio>

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE(struct OrbImplementation;);
VECGEOM_DEVICE_DECLARE_CONV(struct, OrbImplementation);

inline namespace VECGEOM_IMPL_NAMESPACE {

class PlacedOrb;
template <typename T>
struct OrbStruct;
class UnplacedOrb;

/// @brief Implements scalar navigation kernels for `UnplacedOrb`.
struct OrbImplementation {

  using PlacedShape_t    = PlacedOrb;
  using UnplacedStruct_t = OrbStruct<Precision>;
  using UnplacedVolume_t = UnplacedOrb;

  /// @brief Test whether a local point is contained in or on the orb.
  /// @tparam Real_v Floating-point scalar type.
  /// @param orb Orb data containing the radius.
  /// @param point Local point to test.
  /// @param inside Set to true unless the point is outside the outer tolerance band.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void Contains(UnplacedStruct_t const &orb,
                                                                    Vector3D<Real_v> const &point, bool &inside)
  {
    bool unused = false, outside = false;
    GenericKernelForContainsAndInside<Real_v, false>(orb, point, unused, outside);
    inside = !outside;
  }

  /// @brief Classify a local point as inside, outside, or surface.
  /// @tparam Real_v Floating-point scalar type.
  /// @param orb Orb data containing the radius.
  /// @param point Local point to classify.
  /// @param inside Set to `kInside`, `kOutside`, or `kSurface`.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void Inside(UnplacedStruct_t const &orb,
                                                                  Vector3D<Real_v> const &point, Inside_t &inside)
  {

    bool completelyinside = false, completelyoutside = false;
    GenericKernelForContainsAndInside<Real_v, true>(orb, point, completelyinside, completelyoutside);
    inside = EInside::kSurface;
    if (completelyoutside) inside = Inside_t(EInside::kOutside);
    if (completelyinside) inside = Inside_t(EInside::kInside);
  }

  /// @brief Compute strict inside/outside flags for point classification.
  /// @tparam Real_v Floating-point scalar type.
  /// @tparam ForInside When true, also compute the strict-inside flag.
  /// @param orb Orb data containing the radius.
  /// @param localPoint Local point to classify.
  /// @param completelyinside Set when the point is inside the inner tolerance radius.
  /// @param completelyoutside Set when the point is outside the outer tolerance radius.
  ///
  /// @details
  /// This helper compares squared radii because classification only needs an
  /// ordering against radial tolerance limits.
  template <typename Real_v, bool ForInside>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void GenericKernelForContainsAndInside(
      UnplacedStruct_t const &orb, Vector3D<Real_v> const &localPoint, bool &completelyinside, bool &completelyoutside)
  {
    Precision fR = orb.fR;
    Real_v rad2  = localPoint.Mag2();
    Real_v tolR  = fR - Real_v(kTolerance);
    if (ForInside) completelyinside = (rad2 <= tolR * tolR);
    tolR              = fR + Real_v(kTolerance);
    completelyoutside = (rad2 >= tolR * tolR);
    return;
  }

  /// @brief Compute the distance from outside the orb to the first entry.
  /// @tparam Real_v Floating-point scalar type.
  /// @param orb Orb data containing the radius.
  /// @param point Local start point.
  /// @param direction Unit local direction.
  /// @param distance Set to the entry distance, `kInfLength`, or `-1` when the point is clearly inside.
  ///
  /// @details
  /// `stepMax` is intentionally ignored. The setup uses squared radii for the
  /// inside/surface predicates, then passes the precomputed radius squared and
  /// radial projection to the quadratic helper.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void DistanceToIn(UnplacedStruct_t const &orb,
                                                                        Vector3D<Real_v> const &point,
                                                                        Vector3D<Real_v> const &direction,
                                                                        Real_v const & /*stepMax*/, Real_v &distance)
  {
    distance          = kInfLength;
    Real_v rad2       = point.Mag2();
    Real_v tolRadius  = Real_v(orb.fR - kTolerance);
    Real_v tolRadius2 = tolRadius * tolRadius;
    if (rad2 < tolRadius2) {
      distance = Real_v(-1.);
      return;
    }

    Real_v pDotV3D        = point.Dot(direction);
    tolRadius             = Real_v(orb.fR + kTolerance);
    tolRadius2            = tolRadius * tolRadius;
    bool isPointOnSurface = (rad2 <= tolRadius2);
    if (isPointOnSurface) {
      if (pDotV3D < Real_v(-kTolerance)) distance = Real_v(0.);
      return;
    }

    DetectIntersectionAndCalculateDistance<Real_v, true>(orb, rad2, pDotV3D, distance);
  }

  /// @brief Store the supporting-plane result and requested radial normal.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SetExitHit(UnplacedStruct_t const &orb,
                                                                      Vector3D<Real_v> const &point,
                                                                      Vector3D<Real_v> const &direction,
                                                                      Real_v const &distance,
                                                                      SurfaceHitView<Real_v> *hit_info)
  {
    if (hit_info->WantsNormal()) {
      Vector3D<Real_v> hitPoint = point + distance * direction;
      hitPoint *= distance == Real_v(0.) ? Real_v(1.) / Sqrt(hitPoint.Mag2()) : Real_v(1.) / Real_v(orb.fR);
      hit_info->SetNormal(hitPoint);
    }
    hit_info->SetLeftBehind(true);
  }

  /// @brief Execute the exit query with compile-time control of hit recording.
  /// @tparam RecordHitInfo Whether to populate the optional exit information.
  template <typename Real_v, bool RecordHitInfo>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void DistanceToOutImpl(UnplacedStruct_t const &orb,
                                                                             Vector3D<Real_v> const &point,
                                                                             Vector3D<Real_v> const &direction,
                                                                             Real_v &distance,
                                                                             SurfaceHitView<Real_v> *hit_info)
  {
    if constexpr (RecordHitInfo) hit_info->Clear();
    distance = kInfLength;

    Real_v rad2       = point.Mag2();
    Real_v tolRadius  = Real_v(orb.fR + kTolerance);
    Real_v tolRadius2 = tolRadius * tolRadius;
    bool validStart   = true;
    if constexpr (RecordHitInfo) validStart = rad2 < tolRadius2;
    if (rad2 > tolRadius2) {
      distance = Real_v(-1.);
      return;
    }

    Real_v pDotV3D        = point.Dot(direction);
    tolRadius             = Real_v(orb.fR - kTolerance);
    tolRadius2            = tolRadius * tolRadius;
    bool isPointOnSurface = (rad2 >= tolRadius2);
    if (isPointOnSurface) {
      if (pDotV3D >= Real_v(-kTolerance)) {
        distance = Real_v(0.);
        if constexpr (RecordHitInfo) {
          if (validStart) SetExitHit(orb, point, direction, distance, hit_info);
        }
        return;
      }
    }

    const bool crossed = DetectIntersectionAndCalculateDistance<Real_v, false>(orb, rad2, pDotV3D, distance);
    if constexpr (RecordHitInfo) {
      // A tolerated exterior start can still miss the mathematical sphere.
      if (validStart && crossed) SetExitHit(orb, point, direction, distance, hit_info);
    }
  }

  /// @brief Compute the distance from inside the orb to the first exit.
  /// @tparam Real_v Floating-point scalar type.
  /// @param orb Orb data containing the radius.
  /// @param point Local start point.
  /// @param direction Unit local direction.
  /// @param stepMax Ignored distance limit.
  /// @param distance Set to the exit distance, `kInfLength`, or `-1` when the point is clearly outside.
  /// @param hit_info Optional scalar exit-normal destination and supporting-plane flag.
  ///
  /// @details
  /// `stepMax` is intentionally ignored. The setup uses squared radii for the
  /// outside/surface predicates, then passes the precomputed radius squared and
  /// radial projection to the quadratic helper. A regular requested normal
  /// reuses the hit distance and known radius; only a zero-distance surface
  /// start requires explicit normalization. The compile-time split keeps the
  /// ordinary query identical to the pre-existing distance kernel.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void DistanceToOut(UnplacedStruct_t const &orb,
                                                                         Vector3D<Real_v> const &point,
                                                                         Vector3D<Real_v> const &direction,
                                                                         Real_v const & /* stepMax */, Real_v &distance,
                                                                         SurfaceHitView<Real_v> *hit_info = nullptr)
  {
    if (hit_info)
      DistanceToOutImpl<Real_v, true>(orb, point, direction, distance, hit_info);
    else
      DistanceToOutImpl<Real_v, false>(orb, point, direction, distance, nullptr);
  }

  /// @brief Compute safety from an outside point to the orb.
  /// @tparam Real_v Floating-point scalar type.
  /// @param orb Orb data containing the radius.
  /// @param point Local point.
  /// @param safety Set to radial safety, zero in the surface band, or `-1` when clearly inside.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SafetyToIn(UnplacedStruct_t const &orb,
                                                                      Vector3D<Real_v> const &point, Real_v &safety)
  {
    Real_v rad = point.Mag();
    safety     = rad - Real_v(orb.fR);
    if (rad < Real_v(orb.fR - kTolerance)) {
      safety = Real_v(-1.);
      return;
    }

    if (rad > Real_v(orb.fR - kTolerance) && rad < Real_v(orb.fR + kTolerance)) safety = Real_v(0.);
  }

  /// @brief Compute safety from an inside point to the orb boundary.
  /// @tparam Real_v Floating-point scalar type.
  /// @param orb Orb data containing the radius.
  /// @param point Local point.
  /// @param safety Set to radial safety, zero in the surface band, or `-1` when clearly outside.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SafetyToOut(UnplacedStruct_t const &orb,
                                                                       Vector3D<Real_v> const &point, Real_v &safety)
  {
    Real_v rad = point.Mag();
    safety     = Real_v(orb.fR) - rad;

    if (rad > Real_v(orb.fR + kTolerance)) {
      safety = Real_v(-1.);
      return;
    }

    if (rad > Real_v(orb.fR - kTolerance) && rad < Real_v(orb.fR + kTolerance)) safety = Real_v(0.);
  }

  /// @brief Solve the radial ray-sphere intersection.
  /// @tparam Real_v Floating-point scalar type.
  /// @tparam ForDistanceToIn Selects the entry root when true, exit root otherwise.
  /// @param orb Orb data containing the radius.
  /// @param rad2 Squared radius of the start point.
  /// @param pDotV3D Radial projection of the direction.
  /// @param distance Set to the selected intersection distance when a valid root exists.
  /// @return True when a valid root was selected.
  ///
  /// @details
  /// The quadratic is `t^2 + 2 pDotV3D t + (rad2 - r^2) = 0`. Entry uses the
  /// smaller root and requires inward motion; exit uses the larger root.
  template <typename Real_v, bool ForDistanceToIn>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static bool DetectIntersectionAndCalculateDistance(
      UnplacedStruct_t const &orb, Real_v const &rad2, Real_v const &pDotV3D, Real_v &distance)
  {

    Real_v radius = Real_v(orb.fR);
    Real_v c      = rad2 - radius * radius;
    Real_v d2     = (pDotV3D * pDotV3D - c);

    if (ForDistanceToIn) {
      if (d2 < Real_v(0.) || pDotV3D > Real_v(0.)) return false;
      distance = -pDotV3D - Sqrt(d2);
      return true;
    } else {
      if (d2 < Real_v(0.)) return false;
      distance = -pDotV3D + Sqrt(d2);
      return true;
    }
  }

  /// @brief Compute the outward radial normal.
  /// @tparam Real_v Floating-point scalar type.
  /// @param orb Orb data containing the radius.
  /// @param point Local point.
  /// @param valid Set when the point is in the radial surface tolerance band.
  /// @return Unit radial normal at `point`.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static Vector3D<Real_v> NormalKernel(UnplacedStruct_t const &orb,
                                                                                    Vector3D<Real_v> const &point,
                                                                                    bool &valid)
  {
    Real_v rad2             = point.Mag2();
    Real_v invRadius        = Real_v(1.) / Sqrt(rad2);
    Vector3D<Real_v> normal = point * invRadius;

    Real_v tolRMaxO = orb.fR + kTolerance;
    Real_v tolRMaxI = orb.fR - kTolerance;

    valid = ((rad2 <= tolRMaxO * tolRMaxO) && (rad2 >= tolRMaxI * tolRMaxI));
    return normal;
  }
};
} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif // VECGEOM_VOLUMES_KERNEL_orbIMPLEMENTATION_H_
