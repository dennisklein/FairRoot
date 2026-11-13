/*
 * ConeImplementation.h
 *
 *  Created on: May 14, 2014
 *      Author: swenzel
 */

/// History notes:
/// revision + moving to Vectorized Cone Kernels (Raman Sehgal)
/// May-June 2017: revision + moving to new Structure (Raman Sehgal)
/// 20180323 Guilherme Lima  Adapted to new UnplacedVolume factory

#ifndef VECGEOM_VOLUMES_KERNEL_CONEIMPLEMENTATION_H_
#define VECGEOM_VOLUMES_KERNEL_CONEIMPLEMENTATION_H_

#include "VecGeom/base/Vector3D.h"
#include "VecGeom/volumes/kernel/GenericKernels.h"
#include "VecGeom/volumes/kernel/shapetypes/ConeTypes.h"
#include "VecGeom/volumes/ConeStruct.h"
#include "VecGeom/volumes/ConeUtilities.h"
#include "VecGeom/volumes/SurfaceHitView.h"

namespace vecgeom {

VECGEOM_DEVICE_DECLARE_CONV_TEMPLATE_2t(struct, ConeImplementation, typename, typename);

inline namespace VECGEOM_IMPL_NAMESPACE {

template <typename T>
class SPlacedCone;
template <typename T>
class SUnplacedCone;

template <typename coneTypeT, typename ConeStructT = ConeStruct<Precision>>
struct ConeImplementation {

  using UnplacedStruct_t = ConeStructT;
  using UnplacedVolume_t = SUnplacedCone<coneTypeT>;
  using PlacedShape_t    = SPlacedCone<UnplacedVolume_t>;

  enum class ExitSurface : unsigned char { kNone, kZMin, kZMax, kInner, kOuter, kPhiStart, kPhiEnd };

  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SetZExitHit(SurfaceHitView<Real_v> *hit_info, bool upper)
  {
    hit_info->SetLeftBehind(true);
    if (hit_info->WantsNormal()) {
      const Real_v zSign = upper ? Real_v(1.) : Real_v(-1.);
      hit_info->SetNormal(Vector3D<Real_v>(Real_v(0.), Real_v(0.), zSign));
    }
  }

  template <typename Real_v, bool ForInnerSurface>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SetConicalExitHit(UnplacedStruct_t const &cone,
                                                                             Vector3D<Real_v> const &point,
                                                                             Vector3D<Real_v> const &direction,
                                                                             Real_v const &distance,
                                                                             SurfaceHitView<Real_v> *hit_info)
  {
    hit_info->SetLeftBehind(!ForInnerSurface);
    if (hit_info->WantsNormal()) {
      if constexpr (!ForInnerSurface) {
        // The cached secant normalizes the conical slope; only rho needs a square root.
        const Vector3D<Real_v> hit = point + distance * direction;
        const Real_v invRhoSec     = cone.fInvSecRMax / NonZero(hit.Perp());
        hit_info->SetNormal(Vector3D<Real_v>(hit.x() * invRhoSec, hit.y() * invRhoSec, cone.fZNormOuter));
      }
    }
  }

  template <typename Real_v, bool StartSurface>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SetPhiExitHit(UnplacedStruct_t const &cone,
                                                                         SurfaceHitView<Real_v> *hit_info)
  {
    bool const leftBehind = cone.fDPhi <= kPi;
    hit_info->SetLeftBehind(leftBehind);
    if (hit_info->WantsNormal()) {
      if (leftBehind) {
        const Vector3D<Precision> normal = StartSurface ? cone.fPhiWedge.GetNormal1() : cone.fPhiWedge.GetNormal2();
        hit_info->SetNormal(Vector3D<Real_v>(-normal.x(), -normal.y(), Real_v(0.)));
      }
    }
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SetExitHit(UnplacedStruct_t const &cone,
                                                                      Vector3D<Real_v> const &point,
                                                                      Vector3D<Real_v> const &direction,
                                                                      Real_v const &distance, ExitSurface surface,
                                                                      SurfaceHitView<Real_v> *hit_info)
  {
    switch (surface) {
    case ExitSurface::kZMin:
      SetZExitHit(hit_info, false);
      break;
    case ExitSurface::kZMax:
      SetZExitHit(hit_info, true);
      break;
    case ExitSurface::kInner:
      SetConicalExitHit<Real_v, true>(cone, point, direction, distance, hit_info);
      break;
    case ExitSurface::kOuter:
      SetConicalExitHit<Real_v, false>(cone, point, direction, distance, hit_info);
      break;
    case ExitSurface::kPhiStart:
      if constexpr (ConeTypes::NeedsPhiTreatment<coneTypeT>::value != ConeTypes::kNo)
        SetPhiExitHit<Real_v, true>(cone, hit_info);
      break;
    case ExitSurface::kPhiEnd:
      if constexpr (ConeTypes::NeedsPhiTreatment<coneTypeT>::value != ConeTypes::kNo)
        SetPhiExitHit<Real_v, false>(cone, hit_info);
      break;
    case ExitSurface::kNone:
      break;
    }
  }

  /* Check whether a point already known to lie on a z plane is also on the matching ring edge. */
  template <typename Real_v, bool ForInnerSurface, bool ForLowerZ>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static bool IsOnRing(UnplacedStruct_t const &cone,
                                                                    Vector3D<Real_v> const &point)
  {
    Real_v rad2 = point.Perp2();

    if (ForLowerZ) {
      if (ForInnerSurface) {
        return (rad2 <= MakePlusTolerantSquare<true>(cone.fRmin1)) &&
               (rad2 >= MakeMinusTolerantSquare<true>(cone.fRmin1));
      } else {
        return (rad2 <= MakePlusTolerantSquare<true>(cone.fRmax1)) &&
               (rad2 >= MakeMinusTolerantSquare<true>(cone.fRmax1));
      }
    } else {
      if (ForInnerSurface) {
        return (rad2 <= MakePlusTolerantSquare<true>(cone.fRmin2)) &&
               (rad2 >= MakeMinusTolerantSquare<true>(cone.fRmin2));
      } else {
        return (rad2 <= MakePlusTolerantSquare<true>(cone.fRmax2)) &&
               (rad2 >= MakeMinusTolerantSquare<true>(cone.fRmax2));
      }
    }
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static bool IsOnPhiBoundaryRay(Vector3D<Real_v> const &point,
                                                                              Precision alongX, Precision alongY)
  {
    const Real_v planeDistance = (-point.x() * alongY) + (point.y() * alongX);
    return ((point.x() * alongX) + (point.y() * alongY) >= Real_v(0.)) &&
           (Abs(planeDistance) <= Real_v(kConeTolerance));
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static bool IsLeavingOtherBoundaryForPhiRoot(
      UnplacedStruct_t const &cone, Vector3D<Real_v> const &point, Vector3D<Real_v> const &dir, Real_v const &rsq,
      Real_v const &outerRad)
  {
    using namespace ConeUtilities;

    // Suppress a same-boundary phi root only when some other active boundary is
    // already the actual exit direction at the starting edge point.
    bool leavingOtherBoundary = (rsq > MakeMinusTolerantSquare<true>(outerRad, cone.fOuterTolerance)) &&
                                (dir.Dot(GetNormal<Real_v, false>(cone, point)) >= -kHalfTolerance);
    if (ConeTypes::checkRminTreatment<coneTypeT>(cone)) {
      Real_v innerRad     = GetRadiusOfConeAtPoint<Real_v, true>(cone, point.z());
      bool onInnerSurface = rsq < MakePlusTolerantSquare<true>(innerRad, cone.fInnerTolerance);
      leavingOtherBoundary |= onInnerSurface && (dir.Dot(GetNormal<Real_v, true>(cone, point)) >= kHalfTolerance);
    }

    Real_v zSurfaceDistance = Abs(point.z()) - cone.fDz;
    leavingOtherBoundary |= (Abs(zSurfaceDistance) < kConeTolerance) && ((point.z() * dir.z()) > Real_v(0.));
    return leavingOtherBoundary;
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static bool IsOnDegenerateClosingRing(UnplacedStruct_t const &cone,
                                                                                     Vector3D<Real_v> const &point)
  {
    using namespace ConeTypes;
    if (!checkRminTreatment<coneTypeT>(cone)) return false;
    if (Abs(Abs(point.z()) - cone.fDz) > Real_v(kTolerance)) return false;

    const bool onUpperEnd = point.z() >= Real_v(0.);
    const Precision rmin  = onUpperEnd ? cone.fRmin2 : cone.fRmin1;
    const Precision rmax  = onUpperEnd ? cone.fRmax2 : cone.fRmax1;
    if (Abs(rmax - rmin) >= kTolerance) return false;

    const Real_v rsq = point.Perp2();
    return (rsq <= MakePlusTolerantSquare<true>(rmax, cone.fOuterTolerance)) &&
           (rsq >= MakeMinusTolerantSquare<true>(rmax, cone.fOuterTolerance));
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static bool IsEnteringDegenerateClosingRing(
      UnplacedStruct_t const &cone, Vector3D<Real_v> const &point, Vector3D<Real_v> const &dir)
  {
    if (!IsOnDegenerateClosingRing(cone, point)) return false;

    return (dir.Dot(ConeUtilities::GetNormal<Real_v, false>(cone, point)) < -kHalfTolerance) &&
           (dir.Dot(ConeUtilities::GetNormal<Real_v, true>(cone, point)) < -kHalfTolerance);
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static bool RejectDegenerateClosingRingDistanceToInCandidate(
      UnplacedStruct_t const &cone, Vector3D<Real_v> const &point, Vector3D<Real_v> const &dir, Real_v const &distance)
  {
    if (!(distance >= Real_v(0.)) || !(distance < Real_v(kInfLength))) return false;

    const Vector3D<Real_v> hit = point + distance * dir;
    return IsOnDegenerateClosingRing(cone, hit) && !IsEnteringDegenerateClosingRing(cone, hit, dir);
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void Contains(UnplacedStruct_t const &cone,
                                                                    Vector3D<Real_v> const &point, bool &inside)
  {
    bool unused  = false;
    bool outside = false;
    ConeHelpers<Real_v, coneTypeT, UnplacedStruct_t>::template GenericKernelForContainsAndInside<false>(
        cone, point, unused, outside);
    inside = !outside;
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void Inside(UnplacedStruct_t const &cone,
                                                                  Vector3D<Real_v> const &point, Inside_t &inside)
  {
    ConeHelpers<Real_v, coneTypeT, UnplacedStruct_t>::Inside(cone, point, inside);
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void DistanceToInImpl(
      UnplacedStruct_t const &cone, Vector3D<Real_v> const &point, Vector3D<Real_v> const &dir,
      Real_v const & /*stepMax*/, Real_v &distance, evolution::Wedge const *externalWedge)
  {
    using namespace ConeUtilities;
    using namespace ConeTypes;
    const Real_v zero(0.);

    // Reject rays that are already outside and moving away before doing any
    // intersection work.
    distance = kInfLength;

    // outside of Z range and going away?
    Real_v distz         = Abs(point.z()) - cone.fDz; // avoid a division for now
    bool outZAndGoingOut = (distz > kTolerance && (point.z() * dir.z()) >= zero) ||
                           (Abs(distz) < kTolerance && (point.z() * dir.z()) > zero);
    if (outZAndGoingOut) return;

    // outside or *on* outer cone and going away?
    Real_v outerRad              = GetRadiusOfConeAtPoint<Real_v, false>(cone, point.z());
    Real_v rsq                   = point.Perp2();
    bool outsideOuterAndGoingOut = (rsq > MakeMinusTolerantSquare<true>(outerRad, cone.fOuterTolerance)) &&
                                   (dir.Dot(GetNormal<Real_v, false>(cone, point)) >= -kHalfTolerance);
    if (outsideOuterAndGoingOut) return;

    //=== Next, check all dimensions of the cone: for points inside --> return -1
    distance = Real_v(-1.0);

    // For points inside z-range, return -1
    bool inside = distz < -kTolerance;

    inside &= rsq < MakeMinusTolerantSquare<true>(outerRad, cone.fOuterTolerance);

    if (checkRminTreatment<coneTypeT>(cone)) {
      Real_v innerRad = GetRadiusOfConeAtPoint<Real_v, true>(cone, point.z());
      inside &= rsq > MakePlusTolerantSquare<true>(innerRad, cone.fInnerTolerance);
    }
    if constexpr (NeedsPhiTreatment<coneTypeT>::value != kNo) {
      if (checkPhiTreatment<coneTypeT>(cone) && inside) {
        bool completelyinsidephi(false);
        bool completelyoutsidephi(false);
        ClassifyPointInCyclicalSector<Real_v, coneTypeT>(cone, point.x(), point.y(), completelyinsidephi,
                                                         completelyoutsidephi);
        inside &= completelyinsidephi;
      }
    }
    if (externalWedge && inside) inside &= externalWedge->Inside(point) == EInside::kInside;
    if (inside) return;

    // Check the z-plane candidate before the more expensive conical and phi
    // intersections.
    distance = Real_v(kInfLength);

    if (IsOnDegenerateClosingRing(cone, point) && IsEnteringDegenerateClosingRing(cone, point, dir)) {
      distance = zero;
      return;
    }

    if (distz >= -kTolerance && dir.z() != zero) {
      Real_v distToZ = distz / NonZero(Abs(dir.z()));

      Real_v hitx = point.x() + distToZ * dir.x();
      Real_v hity = point.y() + distToZ * dir.y();

      Real_v r2 = (hitx * hitx) + (hity * hity);

      Precision innerZTol       = cone.fTolIz;
      bool isHittingTopPlane    = (point.z() >= innerZTol) && (r2 <= cone.fSqRmax2 + kTolerance);
      bool isHittingBottomPlane = (point.z() <= -innerZTol) && (r2 <= cone.fSqRmax1 + kTolerance);
      bool okz                  = isHittingTopPlane || isHittingBottomPlane;

      if (checkRminTreatment<coneTypeT>(cone)) {
        isHittingTopPlane &= (r2 >= cone.fSqRmin2 - kTolerance);
        isHittingBottomPlane &= (r2 >= cone.fSqRmin1 - kTolerance);
        okz &= (isHittingTopPlane || isHittingBottomPlane);
      }

      if constexpr (NeedsPhiTreatment<coneTypeT>::value != kNo) {
        if (checkPhiTreatment<coneTypeT>(cone) && okz) {
          bool insector(false);
          PointInCyclicalSector<Real_v, coneTypeT, false, true>(cone, hitx, hity, insector, kTolerance);
          okz = insector;
        }
      }
      if (externalWedge && okz) okz = externalWedge->ContainsWithBoundary(Vector3D<Real_v>(hitx, hity, zero));
      if (okz && RejectDegenerateClosingRingDistanceToInCandidate(cone, point, dir, distToZ)) okz = false;
      if (okz) {
        distance = distToZ;
        return;
      }
    }

    Real_v dist_rOuter(kInfLength);
    bool ok_outerCone = ConeHelpers<Real_v, coneTypeT, UnplacedStruct_t>::
        template DetectIntersectionAndCalculateDistanceToConicalSurface<true, false>(cone, point, dir, dist_rOuter,
                                                                                     externalWedge);
    if (ok_outerCone && RejectDegenerateClosingRingDistanceToInCandidate(cone, point, dir, dist_rOuter)) {
      ok_outerCone = false;
    }
    if (ok_outerCone && dist_rOuter < distance) {
      distance = dist_rOuter;
      return;
    }

    Real_v dist_rInner(kInfLength);
    if (checkRminTreatment<coneTypeT>(cone)) {
      bool ok_innerCone = ConeHelpers<Real_v, coneTypeT, UnplacedStruct_t>::
          template DetectIntersectionAndCalculateDistanceToConicalSurface<true, true>(cone, point, dir, dist_rInner,
                                                                                      externalWedge);
      if (ok_innerCone && RejectDegenerateClosingRingDistanceToInCandidate(cone, point, dir, dist_rInner)) {
        ok_innerCone = false;
      }
      if (ok_innerCone && dist_rInner < distance) distance = dist_rInner;
    }

    if constexpr (NeedsPhiTreatment<coneTypeT>::value != kNo) {
      if (checkPhiTreatment<coneTypeT>(cone)) {
        evolution::Wedge const &w = cone.fPhiWedge;

        Real_v dist_phi(kInfLength);
        bool ok_phi(false);
        PhiPlaneTrajectoryIntersection<Real_v, coneTypeT, SectorType<coneTypeT>::value != kOnePi, true>(
            cone.fAlongPhi1x, cone.fAlongPhi1y, w.GetNormal1().x(), w.GetNormal1().y(), cone, point, dir, dist_phi,
            ok_phi);
        if (ok_phi && dist_phi < distance) {
          bool ignoreStartPhiRoot = false;
          if (dist_phi <= Real_v(kConeTolerance)) {
            if (IsOnPhiBoundaryRay(point, cone.fAlongPhi1x, cone.fAlongPhi1y))
              ignoreStartPhiRoot = IsLeavingOtherBoundaryForPhiRoot(cone, point, dir, rsq, outerRad);
          }
          if (!ignoreStartPhiRoot) {
            distance = dist_phi;
            return;
          }
        }

        if (SectorType<coneTypeT>::value != kOnePi) {
          PhiPlaneTrajectoryIntersection<Real_v, coneTypeT, true, true>(cone.fAlongPhi2x, cone.fAlongPhi2y,
                                                                        w.GetNormal2().x(), w.GetNormal2().y(), cone,
                                                                        point, dir, dist_phi, ok_phi);
          if (ok_phi && dist_phi < distance) {
            bool ignoreEndPhiRoot = false;
            if (dist_phi <= Real_v(kConeTolerance)) {
              if (IsOnPhiBoundaryRay(point, cone.fAlongPhi2x, cone.fAlongPhi2y))
                ignoreEndPhiRoot = IsLeavingOtherBoundaryForPhiRoot(cone, point, dir, rsq, outerRad);
            }
            if (!ignoreEndPhiRoot) distance = dist_phi;
          }
        }
      }
    }
    if (externalWedge) {
      const Vector3D<Precision> along1 = externalWedge->GetAlong1();
      const Vector3D<Precision> along2 = externalWedge->GetAlong2();
      Real_v dist_phi1(kInfLength), dist_phi2(kInfLength);
      bool ok_phi1 = false, ok_phi2 = false;
      PhiPlaneTrajectoryIntersection<Real_v, UniversalCone, true, true>(
          along1.x(), along1.y(), externalWedge->GetNormal1().x(), externalWedge->GetNormal1().y(), cone, point, dir,
          dist_phi1, ok_phi1);
      PhiPlaneTrajectoryIntersection<Real_v, UniversalCone, true, true>(
          along2.x(), along2.y(), externalWedge->GetNormal2().x(), externalWedge->GetNormal2().y(), cone, point, dir,
          dist_phi2, ok_phi2);

      const bool usePhi1    = ok_phi1 && (!ok_phi2 || dist_phi1 <= dist_phi2);
      const Real_v dist_phi = usePhi1 ? dist_phi1 : dist_phi2;
      if ((ok_phi1 || ok_phi2) && dist_phi < distance) {
        bool rejectZeroPhiRoot = false;
        if (dist_phi <= Real_v(kConeTolerance)) {
          const bool onPhiBoundary = usePhi1 ? IsOnPhiBoundaryRay(point, along1.x(), along1.y())
                                             : IsOnPhiBoundaryRay(point, along2.x(), along2.y());
          if (onPhiBoundary) rejectZeroPhiRoot = IsLeavingOtherBoundaryForPhiRoot(cone, point, dir, rsq, outerRad);
        }
        if (!rejectZeroPhiRoot) distance = dist_phi;
      }
    }
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void DistanceToIn(UnplacedStruct_t const &cone,
                                                                        Vector3D<Real_v> const &point,
                                                                        Vector3D<Real_v> const &dir,
                                                                        Real_v const &stepMax, Real_v &distance)
  {
    DistanceToInImpl(cone, point, dir, stepMax, distance, nullptr);
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void DistanceToIn(UnplacedStruct_t const &cone,
                                                                        Vector3D<Real_v> const &point,
                                                                        Vector3D<Real_v> const &dir,
                                                                        Real_v const &stepMax, Real_v &distance,
                                                                        evolution::Wedge const *externalWedge)
  {
    DistanceToInImpl(cone, point, dir, stepMax, distance, externalWedge);
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static bool IsIncidentOnZSurface(UnplacedStruct_t const &cone,
                                                                                Vector3D<Real_v> const &point)
  {
    return Abs(Abs(point.z()) - Real_v(cone.fDz)) < Real_v(kConeTolerance);
  }

  template <typename Real_v, bool RecordHitInfo, bool RecordZHit>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static bool DistanceToOutImpl(
      UnplacedStruct_t const &cone, Vector3D<Real_v> const &point, Vector3D<Real_v> const &direction,
      Real_v const & /*stepMax*/, Real_v &distance, SurfaceHitView<Real_v> *hit_info)
  {
    using namespace ConeUtilities;
    using namespace ConeTypes;

    const Real_v zero(0.0);
    if constexpr (RecordHitInfo) hit_info->Clear();
    bool hitZSurface = false;
    [[maybe_unused]] ExitSurface hitSurface{ExitSurface::kNone};

    // Delay the z division until a z-plane candidate is still plausible.
    Real_v distz = Abs(point.z()) - cone.fDz;

    // Reject points that are already outside before looking for an exit.
    distance = Real_v(-1.0);

    bool outside = distz > Real_v(kConeTolerance);

    Real_v rsq      = point.Perp2();
    Real_v outerRad = ConeUtilities::GetRadiusOfConeAtPoint<Real_v, false>(cone, point.z());
    outside |= rsq > MakePlusTolerantSquare<true>(outerRad, cone.fOuterTolerance);
    if (outside) return hitZSurface;

    // rejection of points on surface and exiting
    const bool onDegenerateClosingRing = IsOnDegenerateClosingRing(cone, point);
    const bool enteringDegenerateClosingRing =
        onDegenerateClosingRing && IsEnteringDegenerateClosingRing(cone, point, direction);
    if (onDegenerateClosingRing && !enteringDegenerateClosingRing) {
      distance = zero;
      if constexpr (RecordHitInfo) {
        SetZExitHit(hit_info, point.z() >= zero);
      }
      if constexpr (RecordZHit) hitZSurface = true;
      return hitZSurface;
    }

    bool onSurfaceAndExiting = (rsq > MakeMinusTolerantSquare<true>(outerRad, cone.fOuterTolerance)) &&
                               (direction.Dot(GetNormal<Real_v, false>(cone, point)) > -kHalfTolerance);
    if (!onDegenerateClosingRing && onSurfaceAndExiting) {
      distance = zero;
      if constexpr (RecordHitInfo) SetConicalExitHit<Real_v, false>(cone, point, direction, distance, hit_info);
      if constexpr (RecordZHit) hitZSurface = IsIncidentOnZSurface(cone, point);
      return hitZSurface;
    }

    bool skipRmin = false;
    if (checkRminTreatment<coneTypeT>(cone)) {
      Real_v innerRad = ConeUtilities::GetRadiusOfConeAtPoint<Real_v, true>(cone, point.z());
      if (rsq < MakeMinusTolerantSquare<true>(innerRad, cone.fInnerTolerance)) return hitZSurface;

      // rejection of points on surface and exiting
      bool onInnerSurface = rsq < MakePlusTolerantSquare<true>(innerRad, cone.fInnerTolerance);
      onSurfaceAndExiting = onInnerSurface && (direction.Dot(GetNormal<Real_v, true>(cone, point)) >= kHalfTolerance);
      if (!onDegenerateClosingRing && onSurfaceAndExiting) {
        distance = zero;
        if constexpr (RecordHitInfo) SetConicalExitHit<Real_v, true>(cone, point, direction, distance, hit_info);
        if constexpr (RecordZHit) hitZSurface = IsIncidentOnZSurface(cone, point);
        return hitZSurface;
      }
      skipRmin = onInnerSurface;
    }

    if constexpr (NeedsPhiTreatment<coneTypeT>::value != kNo) {
      if (checkPhiTreatment<coneTypeT>(cone)) {
        bool completelyInsidePhi(false);
        bool completelyOutsidePhi(false);
        ClassifyPointInCyclicalSector<Real_v, coneTypeT>(cone, point.x(), point.y(), completelyInsidePhi,
                                                         completelyOutsidePhi);
        if (completelyOutsidePhi) return hitZSurface;
      }
    }

    bool isGoingUp   = direction.z() > zero;
    bool isGoingDown = direction.z() < zero;
    bool enteringThroughTopZEdge(false);
    bool enteringThroughBottomZEdge(false);
    // Cap-edge entry is relevant only in the existing z-surface band. Both caps
    // use the same conical-side test, so keep it off the ordinary interior path.
    if (Abs(distz) < kTolerance && !onDegenerateClosingRing) {
      const bool leavingTopZ    = isGoingUp && point.z() > zero;
      const bool leavingBottomZ = isGoingDown && point.z() < zero;
      if (leavingTopZ || leavingBottomZ) {
        bool onOuterConeBoundary =
            Abs(SafeDistanceToConicalSurface<Real_v, false>(cone, point)) <= cone.fOuterTolerance;
        bool enteringThroughZEdge =
            onOuterConeBoundary && (direction.Dot(GetNormal<Real_v, false>(cone, point)) < -kHalfTolerance);
        if (checkRminTreatment<coneTypeT>(cone)) {
          bool onInnerConeBoundary =
              Abs(SafeDistanceToConicalSurface<Real_v, true>(cone, point)) <= cone.fInnerTolerance;
          enteringThroughZEdge |=
              onInnerConeBoundary && (direction.Dot(GetNormal<Real_v, true>(cone, point)) < -kHalfTolerance);
        }
        if (!enteringThroughZEdge) {
          distance = distz;
          if constexpr (RecordHitInfo) SetZExitHit(hit_info, leavingTopZ);
          if constexpr (RecordZHit) hitZSurface = true;
          return hitZSurface;
        }
        enteringThroughTopZEdge    = leavingTopZ;
        enteringThroughBottomZEdge = leavingBottomZ;
      }
    }

    // Check the z-plane candidate first; it can terminate the query before the
    // conical and phi-boundary work.
    distance = Real_v(kInfLength);

    Precision fDz  = cone.fDz;
    Real_v dirZInv = Real_v(1.) / NonZero(direction.z());
    if (isGoingUp && !enteringThroughTopZEdge) {
      distance = (fDz - point.z()) * dirZInv;
      if constexpr (RecordHitInfo) hitSurface = ExitSurface::kZMax;
      if constexpr (RecordZHit) hitZSurface = true;
    }
    if (isGoingDown && !enteringThroughBottomZEdge) {
      distance = (-fDz - point.z()) * dirZInv;
      if constexpr (RecordHitInfo) hitSurface = ExitSurface::kZMin;
      if constexpr (RecordZHit) hitZSurface = true;
    }

    Real_v dist_rOuter(kInfLength);
    bool ok_outerCone = ConeHelpers<Real_v, coneTypeT, UnplacedStruct_t>::
        template DetectIntersectionAndCalculateDistanceToConicalSurface<false, false>(cone, point, direction,
                                                                                      dist_rOuter);
    if (ok_outerCone && dist_rOuter < distance) {
      distance = dist_rOuter;
      if constexpr (RecordHitInfo) hitSurface = ExitSurface::kOuter;
      if constexpr (RecordZHit) hitZSurface = false;
    }

    Real_v dist_rInner(kInfLength);
    if (checkRminTreatment<coneTypeT>(cone) && !skipRmin) {
      bool ok_innerCone = ConeHelpers<Real_v, coneTypeT, UnplacedStruct_t>::
          template DetectIntersectionAndCalculateDistanceToConicalSurface<false, true>(cone, point, direction,
                                                                                       dist_rInner);
      if (ok_innerCone && dist_rInner < distance) {
        distance = dist_rInner;
        if constexpr (RecordHitInfo) hitSurface = ExitSurface::kInner;
        if constexpr (RecordZHit) hitZSurface = false;
      }
    }

    if constexpr (NeedsPhiTreatment<coneTypeT>::value != kNo) {
      if (checkPhiTreatment<coneTypeT>(cone)) {

        bool isOnStartPhi        = ConeUtilities::IsOnStartPhi<Real_v>(cone, point);
        bool isOnEndPhi          = ConeUtilities::IsOnEndPhi<Real_v>(cone, point);
        Vector3D<Real_v> normal1 = cone.fPhiWedge.GetNormal1();
        Vector3D<Real_v> normal2 = cone.fPhiWedge.GetNormal2();
        bool leavingPhiSurface =
            (isOnStartPhi && direction.Dot(-normal1) > zero) || (isOnEndPhi && direction.Dot(-normal2) > zero);
        if (leavingPhiSurface) {
          distance = zero;
          if constexpr (RecordHitInfo) {
            if (isOnStartPhi)
              SetPhiExitHit<Real_v, true>(cone, hit_info);
            else
              SetPhiExitHit<Real_v, false>(cone, hit_info);
          }
          if constexpr (RecordZHit) hitZSurface = IsIncidentOnZSurface(cone, point);
          return hitZSurface;
        }

        Real_v dist_phi(kInfLength);
        bool ok_phi               = false;
        evolution::Wedge const &w = cone.fPhiWedge;
        PhiPlaneTrajectoryIntersection<Real_v, coneTypeT, SectorType<coneTypeT>::value != kOnePi, false>(
            cone.fAlongPhi1x, cone.fAlongPhi1y, w.GetNormal1().x(), w.GetNormal1().y(), cone, point, direction,
            dist_phi, ok_phi);
        if (ok_phi && dist_phi < distance) {
          distance = dist_phi;
          if constexpr (RecordHitInfo) hitSurface = ExitSurface::kPhiStart;
          if constexpr (RecordZHit) hitZSurface = false;
        }

        if (SectorType<coneTypeT>::value != kOnePi) {
          ConeUtilities::PhiPlaneTrajectoryIntersection<Real_v, coneTypeT, true, false>(
              cone.fAlongPhi2x, cone.fAlongPhi2y, w.GetNormal2().x(), w.GetNormal2().y(), cone, point, direction,
              dist_phi, ok_phi);
          if (ok_phi && dist_phi < distance) {
            distance = dist_phi;
            if constexpr (RecordHitInfo) hitSurface = ExitSurface::kPhiEnd;
            if constexpr (RecordZHit) hitZSurface = false;
          }
        }
      }
    }
    if (distance < zero && Abs(distance) < kTolerance) distance = zero;
    if constexpr (RecordHitInfo) {
      if (hitSurface != ExitSurface::kNone) SetExitHit(cone, point, direction, distance, hitSurface, hit_info);
    }
    if constexpr (RecordZHit) {
      if (!hitZSurface) hitZSurface = IsIncidentOnZSurface(cone, point + distance * direction);
    }
    return hitZSurface;
  }

  /**
   * @brief Compute the exit distance and optionally return a supporting normal.
   * @param cone Cone geometry.
   * @param point Starting point in the cone frame.
   * @param direction Normalized query direction.
   * @param stepMax Caller-provided limit; unused by the exit kernel.
   * @param[out] distance Exit distance, or `-1` for a rejected outside point.
   * @param[out] hit_info Optional supporting-plane flag and caller-owned normal storage.
   * @details The selected boundary qualifies the output, including tolerated
   * surface starts. Non-supporting exits leave normal storage untouched.
   */
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void DistanceToOut(UnplacedStruct_t const &cone,
                                                                         Vector3D<Real_v> const &point,
                                                                         Vector3D<Real_v> const &direction,
                                                                         Real_v const &stepMax, Real_v &distance,
                                                                         SurfaceHitView<Real_v> *hit_info = nullptr)
  {
    if (hit_info)
      DistanceToOutImpl<Real_v, true, false>(cone, point, direction, stepMax, distance, hit_info);
    else
      DistanceToOutImpl<Real_v, false, false>(cone, point, direction, stepMax, distance, nullptr);
  }

  /**
   * @brief Compute an exit distance and report whether the hit is incident on a z cap.
   * @details This is a narrow internal protocol for section-based composites such as Polycone. It deliberately
   * avoids producing general surface topology when only z-section hand-off is needed.
   * @param cone Cone section data.
   * @param point Query point in the cone frame.
   * @param direction Normalized query direction in the cone frame.
   * @param stepMax Caller-provided distance limit.
   * @param distance Output exit distance.
   * @return True when the selected exit is incident on either cone z cap.
   */
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static bool DistanceToOutWithZHit(UnplacedStruct_t const &cone,
                                                                                 Vector3D<Real_v> const &point,
                                                                                 Vector3D<Real_v> const &direction,
                                                                                 Real_v const &stepMax,
                                                                                 Real_v &distance)
  {
    return DistanceToOutImpl<Real_v, false, true>(cone, point, direction, stepMax, distance, nullptr);
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SafetyToIn(UnplacedStruct_t const &cone,
                                                                      Vector3D<Real_v> const &point, Real_v &safety)
  {
    using namespace ConeUtilities;
    using namespace ConeTypes;

    safety        = Real_v(-1.0);
    Precision fDz = cone.fDz;
    Real_v distz  = Abs(point.z()) - fDz;

    // For points inside z-range, return -1
    bool inside = distz < -kConeTolerance;

    // This logic to check if the point is inside is far better than
    // using GenericKernel and will improve performance.
    Real_v outerRad = GetRadiusOfConeAtPoint<Real_v, false>(cone, point.z());
    Real_v rsq      = point.Perp2();
    inside &= rsq < MakeMinusTolerantSquare<true>(outerRad, cone.fOuterTolerance);

    Real_v innerRad(0.);
    if (checkRminTreatment<coneTypeT>(cone)) {
      innerRad = GetRadiusOfConeAtPoint<Real_v, true>(cone, point.z());
      inside &= rsq > MakePlusTolerantSquare<true>(innerRad, cone.fInnerTolerance);
    }
    Real_v safetyPhi(-kInfLength);
    if constexpr (NeedsPhiTreatment<coneTypeT>::value != kNo) {
      if (checkPhiTreatment<coneTypeT>(cone) && inside) {
        safetyPhi = cone.fPhiWedge.SafetyToIn(point);
        inside    = safetyPhi < Real_v(0.);
      }
    }
    if (inside) return;

    // Once it is checked that the point is inside or not, safety can be set to 0.
    // This will serve the case that the point is on the surface. So no need to check
    // that the point is really on surface.
    safety = Real_v(0.);

    // Now if the point is neither inside nor on surface, then it should be outside
    // and the safety should be set to some finite value, which is done by below logic

    Real_v safeZ = distz;
    // Reuse the squared radius already computed for the inside test; the conical safety helper needs rho only here.
    Real_v rho                  = Sqrt(rsq);
    Real_v safeDistOuterSurface = -SafeDistanceToConicalSurface<Real_v, false>(cone, point, rho);
    // Reuse the classifier's squared band exactly, including its endpoints.
    if (rsq >= MakeMinusTolerantSquare<true>(outerRad, cone.fOuterTolerance) &&
        rsq <= MakePlusTolerantSquare<true>(outerRad, cone.fOuterTolerance))
      safeDistOuterSurface = Real_v(0.);

    Real_v safeDistInnerSurface(-kInfLength);
    if (checkRminTreatment<coneTypeT>(cone)) {
      safeDistInnerSurface = -SafeDistanceToConicalSurface<Real_v, true>(cone, point, rho);
      if (rsq >= MakeMinusTolerantSquare<true>(innerRad, cone.fInnerTolerance) &&
          rsq <= MakePlusTolerantSquare<true>(innerRad, cone.fInnerTolerance))
        safeDistInnerSurface = Real_v(0.);
    }

    safety = Max(safeZ, Max(safeDistOuterSurface, safeDistInnerSurface));
    if constexpr (NeedsPhiTreatment<coneTypeT>::value != kNo) {
      if (checkPhiTreatment<coneTypeT>(cone)) safety = Max(safetyPhi, safety);
    }

    if (Abs(safety) < kTolerance) safety = Real_v(0.);
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SafetyToOut(UnplacedStruct_t const &cone,
                                                                       Vector3D<Real_v> const &point, Real_v &safety)
  {

    using namespace ConeUtilities;
    using namespace ConeTypes;

    safety = Real_v(-1.0);

    Real_v distz = Abs(point.z()) - cone.fDz;
    Real_v rsq   = point.Perp2();

    // This logic to check if the point is outside is far better then
    // using GenericKernel and will improve performance.
    bool outside = distz > Real_v(kConeTolerance);

    Real_v outerRad = GetRadiusOfConeAtPoint<Real_v, false>(cone, point.z());
    outside |= rsq > MakePlusTolerantSquare<true>(outerRad, cone.fOuterTolerance);

    Real_v innerRad(0.);
    if (checkRminTreatment<coneTypeT>(cone)) {
      innerRad = GetRadiusOfConeAtPoint<Real_v, true>(cone, point.z());
      outside |= rsq < MakeMinusTolerantSquare<true>(innerRad, cone.fInnerTolerance);
    }

    Real_v safetyPhi(kInfLength);
    if constexpr (NeedsPhiTreatment<coneTypeT>::value != kNo) {
      if (checkPhiTreatment<coneTypeT>(cone) && !outside) {
        safetyPhi = cone.fPhiWedge.SafetyToOut(point);
        outside |= safetyPhi < Real_v(0.);
      }
    }
    if (outside) return;

    // Once it is checked that the point is inside or not, safety can be set to 0.
    // This will serve the case that the point is on the surface. So no need to check
    // that the point is really on surface.
    safety = Real_v(0.);

    // Now if the point is neither outside nor on surface, then it should be inside
    // and the safety should be set to some finite value, which is done by below logic

    Precision fDz = cone.fDz;
    Real_v safeZ  = fDz - Abs(point.z());

    Real_v safeDistOuterSurface = SafeDistanceToConicalSurface<Real_v, false>(cone, point);
    // Reuse the classifier's squared band exactly, including its endpoints.
    if (rsq >= MakeMinusTolerantSquare<true>(outerRad, cone.fOuterTolerance) &&
        rsq <= MakePlusTolerantSquare<true>(outerRad, cone.fOuterTolerance))
      safeDistOuterSurface = Real_v(0.);
    Real_v safeDistInnerSurface(kInfLength);
    if (checkRminTreatment<coneTypeT>(cone)) {
      safeDistInnerSurface = SafeDistanceToConicalSurface<Real_v, true>(cone, point);
      if (rsq >= MakeMinusTolerantSquare<true>(innerRad, cone.fInnerTolerance) &&
          rsq <= MakePlusTolerantSquare<true>(innerRad, cone.fInnerTolerance))
        safeDistInnerSurface = Real_v(0.);
    }

    safety = Min(safeZ, Min(safeDistOuterSurface, safeDistInnerSurface));
    if constexpr (NeedsPhiTreatment<coneTypeT>::value != kNo) {
      if (checkPhiTreatment<coneTypeT>(cone)) safety = Min(safetyPhi, safety);
    }
    if (Abs(safety) < kTolerance) safety = Real_v(0.);
  }

  template <typename Real_v, bool ForInnerSurface>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static Real_v SafeDistanceToConicalSurface(UnplacedStruct_t const &cone,
                                                                                          Vector3D<Real_v> const &point,
                                                                                          Real_v const &rho)
  {

    typedef Real_v Float_t;
    if (ForInnerSurface) {
      Float_t pRMin = cone.fTanRMin * point.z() + (cone.fRmin1 + cone.fRmin2) * Float_t(0.5); // cone.fRminAv;
      return (rho - pRMin) * cone.fInvSecRMin;
    } else {
      Float_t pRMax = cone.fTanRMax * point.z() + (cone.fRmax1 + cone.fRmax2) * Float_t(0.5); // cone.fRmaxAv;
      return (pRMax - rho) * cone.fInvSecRMax;
    }
  }

  template <typename Real_v, bool ForInnerSurface>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static Real_v SafeDistanceToConicalSurface(UnplacedStruct_t const &cone,
                                                                                          Vector3D<Real_v> const &point)
  {
    return SafeDistanceToConicalSurface<Real_v, ForInnerSurface>(cone, point, point.Perp());
  }
};
} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif
