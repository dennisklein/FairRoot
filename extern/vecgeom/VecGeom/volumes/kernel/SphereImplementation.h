
/// @file SphereImplementation.h
/// @brief Navigation kernels for spherical shells with optional phi and theta cuts.
/// @author Raman Sehgal (raman.sehgal@cern.ch)

#ifndef VECGEOM_VOLUMES_KERNEL_SPHEREIMPLEMENTATION_H_
#define VECGEOM_VOLUMES_KERNEL_SPHEREIMPLEMENTATION_H_

#include "VecGeom/base/Vector3D.h"
#include "VecGeom/volumes/SphereStruct.h"
#include "VecGeom/volumes/SurfaceHitView.h"
#include "VecGeom/volumes/kernel/GenericKernels.h"
#include <VecCore/VecCore>
#include "VecGeom/volumes/kernel/OrbImplementation.h"
#include "VecGeom/volumes/SphereUtilities.h"

#include <cstdio>

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE(struct SphereImplementation;);
VECGEOM_DEVICE_DECLARE_CONV(struct, SphereImplementation);

inline namespace VECGEOM_IMPL_NAMESPACE {

class PlacedSphere;
template <typename T>
struct SphereStruct;
class UnplacedSphere;

/// @brief Implementation for full, hollow, and angular-cut spheres.
///
/// @details A sphere is modeled as a radial interval clipped by optional phi
/// and theta boundaries. Full variants skip the corresponding angular checks.
struct SphereImplementation {

  using PlacedShape_t    = PlacedSphere;
  using UnplacedStruct_t = SphereStruct<Precision>;
  using UnplacedVolume_t = UnplacedSphere;

  /// Boundary selected by the existing DistanceToOut minimization.
  enum class ExitSurface : unsigned char { kNone, kOuter, kInner, kPhiStart, kPhiEnd, kThetaStart, kThetaEnd };

  /// @brief Compute one phi-plane exit, retaining the wedge's surface and motion tolerances.
  /// @return Forward distance, zero for an outward surface start, or `kInfLength` for no exit.
  template <typename Real_v, bool Start>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static Real_v PhiDistanceToOut(UnplacedStruct_t const &sphere,
                                                                              Vector3D<Real_v> const &point,
                                                                              Vector3D<Real_v> const &direction)
  {
    Vector3D<Precision> const normal = Start ? sphere.fPhiWedge.GetNormal1() : sphere.fPhiWedge.GetNormal2();
    const Real_v motion              = direction.x() * normal.x() + direction.y() * normal.y();
    if (motion >= Real_v(-1.e-12)) return Real_v(kInfLength);

    const Real_v projection         = point.x() * normal.x() + point.y() * normal.y();
    Vector3D<Precision> const along = Start ? sphere.fPhiWedge.GetAlong1() : sphere.fPhiWedge.GetAlong2();
    if (Abs(projection) < Real_v(kTolerance) && point.x() * along.x() + point.y() * along.y() >= Real_v(0.))
      return Real_v(0.);

    Real_v distance = -projection / motion;
    if (distance > Real_v(-kTolerance) && distance < Real_v(0.)) distance = Real_v(0.);
    return distance >= Real_v(0.) ? distance : Real_v(kInfLength);
  }

  /// @brief Compute both phi-plane exit candidates; inactive candidates remain `kInfLength`.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void PhiDistancesToOut(UnplacedStruct_t const &sphere,
                                                                             Vector3D<Real_v> const &point,
                                                                             Vector3D<Real_v> const &direction,
                                                                             Real_v &startDistance, Real_v &endDistance)
  {
    startDistance = PhiDistanceToOut<Real_v, true>(sphere, point, direction);
    endDistance   = PhiDistanceToOut<Real_v, false>(sphere, point, direction);
  }

  /// @brief Compute the exit-directed root of one theta boundary.
  /// @details Sine/cosine scaling keeps the quadratic finite near the equator.
  /// The caller shares the transverse point, direction, and mixed products between both boundaries.
  /// `intersects` excludes the opposite cone nappe and tangencies without an apex crossing.
  template <typename Real_v, bool Start>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void ThetaDistanceToOut(
      UnplacedStruct_t const &sphere, Vector3D<Real_v> const &point, Vector3D<Real_v> const &direction,
      Real_v const &rho2, Real_v const &pDotV2d, Real_v const &dirRho2, Real_v &distance, bool &intersects)
  {
    const Precision theta = Start ? sphere.fSTheta : sphere.eTheta;
    distance              = Real_v(kInfLength);
    intersects            = false;
    if ((Start && theta == 0.) || (!Start && theta == kPi)) return;

    if (Abs(theta - kHalfPi) <= 0.5 * sphere.kAngTolerance) {
      const bool movingOut = Start ? direction.z() > Real_v(0.) : direction.z() < Real_v(0.);
      if (!movingOut) return;
      distance = -point.z() / direction.z();
      if (distance > Real_v(-kTolerance) && distance < Real_v(0.)) distance = Real_v(0.);
      intersects = distance >= Real_v(0.);
      return;
    }

    const Real_v sinTheta = Start ? Real_v(sphere.sinSTheta) : Real_v(sphere.sinETheta);
    const Real_v cosTheta = Start ? Real_v(sphere.cosSTheta) : Real_v(sphere.cosETheta);
    const Real_v sin2     = sinTheta * sinTheta;
    const Real_v cos2     = cosTheta * cosTheta;
    const Real_v a        = cos2 * dirRho2 - sin2 * direction.z() * direction.z();
    const Real_v b        = cos2 * pDotV2d - sin2 * point.z() * direction.z();
    const Real_v c        = cos2 * rho2 - sin2 * point.z() * point.z();
    Real_v discriminant   = b * b - a * c;
    if (discriminant < Real_v(0.) && Abs(discriminant) < Real_v(kHalfTolerance)) discriminant = Real_v(0.);
    if (discriminant < Real_v(0.)) return;

    const Real_v root = Sqrt(discriminant);
    // The derivative at the two roots is +/-2*sqrt(discriminant); its exit sign selects a single root.
    const bool negativeRoot = (Start && theta < kHalfPi) || (!Start && theta > kHalfPi);
    if (negativeRoot) {
      if (b >= Real_v(0.) && a != Real_v(0.))
        distance = (-b - root) / a;
      else if (b < Real_v(0.))
        distance = c / NonZero(-b + root);
    } else {
      if (b > Real_v(0.))
        distance = c / NonZero(-b - root);
      else if (a != Real_v(0.))
        distance = (-b + root) / a;
    }

    if (distance > Real_v(-kTolerance) && distance < Real_v(0.)) distance = Real_v(0.);
    if (distance < Real_v(0.) || distance >= Real_v(kInfLength)) {
      distance = Real_v(kInfLength);
      return;
    }
    const Real_v hitZ    = point.z() + distance * direction.z();
    bool crossesBoundary = discriminant > Real_v(0.);
    if (!crossesBoundary) {
      const Real_v hitX   = point.x() + distance * direction.x();
      const Real_v hitY   = point.y() + distance * direction.y();
      const bool hitsApex = hitX * hitX + hitY * hitY + hitZ * hitZ <= Real_v(kTolerance * kTolerance);
      if (hitsApex)
        crossesBoundary = Start ? direction.z() > Real_v(sphere.cosSTheta) : direction.z() < Real_v(sphere.cosETheta);
    }
    intersects = crossesBoundary && (theta < kHalfPi ? hitZ > Real_v(-kHalfTolerance) : hitZ < Real_v(kHalfTolerance));
  }

  /// @brief Compute both theta exits, sharing the transverse quadratic coefficients.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void ThetaDistancesToOut(
      UnplacedStruct_t const &sphere, Vector3D<Real_v> const &point, Vector3D<Real_v> const &direction,
      Real_v &startDistance, Real_v &endDistance, bool &startIntersects, bool &endIntersects)
  {
    const Real_v rho2    = point.Perp2();
    const Real_v pDotV2d = point.x() * direction.x() + point.y() * direction.y();
    const Real_v dirRho2 = direction.Perp2();
    ThetaDistanceToOut<Real_v, true>(sphere, point, direction, rho2, pDotV2d, dirRho2, startDistance, startIntersects);
    ThetaDistanceToOut<Real_v, false>(sphere, point, direction, rho2, pDotV2d, dirRho2, endDistance, endIntersects);
  }

  /// Return the boundary responsible for an immediate surface exit.
  /// @param sphere Sphere parameters and cached angular state.
  /// @param point Local point classified on the surface.
  /// @param direction Normalized local exit direction.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static ExitSurface SelectSurfaceExit(UnplacedStruct_t const &sphere,
                                                                                    Vector3D<Real_v> const &point,
                                                                                    Vector3D<Real_v> const &direction)
  {
    const Real_v rad2 = point.Mag2();
    if (SphereUtilities::IsPointOnRadialSurfaceAndMovingOut<Real_v, false, true>(sphere, point, direction) &&
        sphere.fPhiWedge.Contains<Real_v>(point) && sphere.fThetaCone.Contains<Real_v>(point))
      return ExitSurface::kOuter;
    if (sphere.fRmin &&
        SphereUtilities::IsPointOnRadialSurfaceAndMovingOut<Real_v, true, true>(sphere, point, direction) &&
        sphere.fPhiWedge.Contains<Real_v>(point) && sphere.fThetaCone.Contains<Real_v>(point))
      return ExitSurface::kInner;
    if (!sphere.fFullThetaSphere && rad2 >= sphere.fRmin * sphere.fRmin && rad2 <= sphere.fRmax * sphere.fRmax &&
        sphere.fPhiWedge.Contains<Real_v>(point)) {
      if (SphereUtilities::IsPointOnThetaSurfaceAndMovingOut<Real_v, true, true>(sphere, point, direction))
        return ExitSurface::kThetaStart;
      if (SphereUtilities::IsPointOnThetaSurfaceAndMovingOut<Real_v, false, true>(sphere, point, direction))
        return ExitSurface::kThetaEnd;
    }
    if (!sphere.fFullPhiSphere && rad2 >= sphere.fRmin * sphere.fRmin && rad2 <= sphere.fRmax * sphere.fRmax &&
        sphere.fThetaCone.Contains<Real_v>(point)) {
      if (sphere.fPhiWedge.IsPointOnSurfaceAndMovingOut<Real_v, true, true>(point, direction))
        return ExitSurface::kPhiStart;
      if (sphere.fPhiWedge.IsPointOnSurfaceAndMovingOut<Real_v, false, true>(point, direction))
        return ExitSurface::kPhiEnd;
    }
    return ExitSurface::kNone;
  }

  /// Convert the recorded winning boundary into the optional exit-hit outputs.
  /// @details A normal is returned only when the selected boundary is a supporting
  /// surface, so the whole solid lies behind its outward tangent plane.
  /// @param sphere Sphere parameters and cached angular state.
  /// @param point Local starting point.
  /// @param direction Normalized local direction.
  /// @param distance Computed exit distance.
  /// @param surface Boundary selected by the distance minimization.
  /// @param hitInfo Optional-output view to fill.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SetExitHit(UnplacedStruct_t const &sphere,
                                                                      Vector3D<Real_v> const &point,
                                                                      Vector3D<Real_v> const &direction,
                                                                      Real_v const &distance, ExitSurface surface,
                                                                      SurfaceHitView<Real_v> *hitInfo)
  {
    bool leftBehind = false;
    switch (surface) {
    case ExitSurface::kOuter:
      leftBehind = true;
      if (hitInfo->WantsNormal()) {
        Vector3D<Real_v> normal = point + distance * direction;
        normal *= distance == Real_v(0.) ? Real_v(1.) / Sqrt(normal.Mag2()) : Real_v(1.) / sphere.fRmax;
        hitInfo->SetNormal(normal);
      }
      break;
    case ExitSurface::kPhiStart:
    case ExitSurface::kPhiEnd:
      leftBehind = sphere.fDPhi <= kPi;
      if (leftBehind && hitInfo->WantsNormal()) {
        const Vector3D<Precision> normal =
            surface == ExitSurface::kPhiStart ? sphere.fPhiWedge.GetNormal1() : sphere.fPhiWedge.GetNormal2();
        hitInfo->SetNormal(Vector3D<Real_v>(-normal.x(), -normal.y(), Real_v(0.)));
      }
      break;
    case ExitSurface::kThetaStart:
    case ExitSurface::kThetaEnd: {
      const bool start = surface == ExitSurface::kThetaStart;
      leftBehind       = start ? sphere.fSTheta >= kHalfPi : sphere.eTheta <= kHalfPi;
      if (leftBehind && hitInfo->WantsNormal()) {
        const Precision theta = start ? sphere.fSTheta : sphere.eTheta;
        Vector3D<Real_v> normal;
        if (Abs(theta - kHalfPi) <= kTolerance) {
          normal.Set(Real_v(0.), Real_v(0.), start ? Real_v(1.) : Real_v(-1.));
        } else {
          Vector3D<Real_v> reference = point + distance * direction;
          Real_v transverse2         = reference.Perp2();
          // Every angular boundary meets at the origin. Use the incident radial direction there because the apex
          // itself has no azimuth from which to construct a cone normal.
          if (transverse2 <= Real_v(kTolerance * kTolerance)) {
            reference   = point.Perp2() > Real_v(kTolerance * kTolerance) ? point : -direction;
            transverse2 = reference.Perp2();
          }
          const Real_v invTransverse = transverse2 > Real_v(0.) ? Real_v(1.) / Sqrt(transverse2) : Real_v(0.);
          const Real_v sinTheta      = start ? sphere.sinSTheta : sphere.sinETheta;
          const Real_v cosTheta      = start ? sphere.cosSTheta : sphere.cosETheta;
          normal.Set(cosTheta * reference.x() * invTransverse, cosTheta * reference.y() * invTransverse, -sinTheta);
          if (transverse2 <= Real_v(0.)) normal.x() = cosTheta;
          if (start) normal = -normal;
        }
        hitInfo->SetNormal(normal);
      }
      break;
    }
    case ExitSurface::kInner:
    case ExitSurface::kNone:
      break;
    }
    hitInfo->SetLeftBehind(leftBehind);
  }

  /// @brief Return whether a point is inside or on the sphere.
  /// @param sphere Sphere parameters and cached angular state.
  /// @param point Local point to test.
  /// @param inside Set to true for inside or surface points.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void Contains(UnplacedStruct_t const &sphere,
                                                                    Vector3D<Real_v> const &point, bool &inside)
  {
    bool unused = false, outside = false;
    GenericKernelForContainsAndInside<Real_v, false>(sphere, point, unused, outside);
    inside = !outside;
  }

  /// @brief Classify a point as inside, outside, or on the surface.
  /// @param sphere Sphere parameters and cached angular state.
  /// @param point Local point to classify.
  /// @param inside Set to the corresponding `EInside` value.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void Inside(UnplacedStruct_t const &sphere,
                                                                  Vector3D<Real_v> const &point, Inside_t &inside)
  {
    bool completelyinside = false, completelyoutside = false;
    GenericKernelForContainsAndInside<Real_v, true, true>(sphere, point, completelyinside, completelyoutside);
    inside = EInside::kSurface;
    if (completelyoutside) inside = Inside_t(EInside::kOutside);
    if (completelyinside) inside = Inside_t(EInside::kInside);
  }

  /// @brief Evaluate complete-inside and complete-outside predicates.
  /// @details Combines the radial interval with active phi and theta clips.
  /// Points in tolerance regions can leave both outputs false.
  /// @tparam ForInside Enables complete-inside updates used by `Inside`.
  /// @tparam EarlyReturn Skips angular work after definite rejection in standalone classification.
  /// @param sphere Sphere parameters and cached angular state.
  /// @param localPoint Local point to classify.
  /// @param completelyinside Set when all active predicates are safely inside.
  /// @param completelyoutside Set when any active predicate is safely outside.
  template <typename Real_v, bool ForInside, bool EarlyReturn = false>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void GenericKernelForContainsAndInside(
      UnplacedStruct_t const &sphere, Vector3D<Real_v> const &localPoint, bool &completelyinside,
      bool &completelyoutside)
  {
    Real_v rad2 = localPoint.Mag2();

    // Check radial surfaces
    // Radial check for GenericKernel Start
    if (sphere.fRmin)
      completelyinside =
          rad2 <= MakeMinusTolerantSquare<true>(sphere.fRmax) && rad2 >= MakePlusTolerantSquare<true>(sphere.fRmin);
    else
      completelyinside = rad2 <= MakeMinusTolerantSquare<true>(sphere.fRmax);

    if (sphere.fRmin)
      completelyoutside =
          rad2 >= MakePlusTolerantSquare<true>(sphere.fRmax) || rad2 <= MakeMinusTolerantSquare<true>(sphere.fRmin);
    else
      completelyoutside = rad2 >= MakePlusTolerantSquare<true>(sphere.fRmax);

    // Phi boundaries  : Do not check if it has no phi boundary!
    if (!sphere.fFullPhiSphere) {
      // Only Inside opts in; overlapping rounded tolerance bounds still need angular checks.
      if constexpr (EarlyReturn) {
        if (completelyoutside && !completelyinside) return;
      }

      bool completelyoutsidephi = false;
      bool completelyinsidephi  = false;
      sphere.fPhiWedge.GenericKernelForContainsAndInside<Real_v, ForInside>(localPoint, completelyinsidephi,
                                                                            completelyoutsidephi);
      completelyoutside |= completelyoutsidephi;

      if (ForInside) completelyinside &= completelyinsidephi;
    }
    // Phi Check for GenericKernel Over

    // Theta bondaries
    if (!sphere.fFullThetaSphere) {
      if constexpr (EarlyReturn) {
        if (completelyoutside && !completelyinside) return;
      }

      bool completelyoutsidetheta = false;
      bool completelyinsidetheta  = false;
      sphere.fThetaCone.GenericKernelForContainsAndInside<Real_v, ForInside>(localPoint, completelyinsidetheta,
                                                                             completelyoutsidetheta);
      completelyoutside |= completelyoutsidetheta;

      if (ForInside) completelyinside &= completelyinsidetheta;
    }
    return;
  }

  /// @brief Compute the first valid entry distance.
  /// @details Radial shell roots, phi-plane roots, and theta-cone roots are
  /// filtered against the active cuts. Theta roots must cross into material;
  /// grazing roots are ignored except for the explicit apex-entry case.
  /// @param sphere Sphere parameters and cached angular state.
  /// @param point Local starting point.
  /// @param direction Normalized local direction.
  /// @param distance Set to the entry distance, `-1`, `0`, or `kInfLength`.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void DistanceToIn(UnplacedStruct_t const &sphere,
                                                                        Vector3D<Real_v> const &point,
                                                                        Vector3D<Real_v> const &direction,
                                                                        Real_v const & /* stepMax */, Real_v &distance)
  {
    distance  = kInfLength;
    bool done = false;

    bool fullPhiSphere   = sphere.fFullPhiSphere;
    bool fullThetaSphere = sphere.fFullThetaSphere;

    // General Precalcs
    Real_v rad2    = point.Mag2();
    Real_v pDotV3d = point.Dot(direction);

    Real_v c = rad2 - sphere.fRmax * sphere.fRmax;

    // Beyond the outer radial tolerance, neither an inside nor a surface start is possible.
    if (rad2 <= MakePlusTolerantSquare<true>(sphere.fRmax)) {
      if (SphereUtilities::IsCompletelyInside<Real_v>(sphere, point)) {
        distance = Real_v(-1.0);
        return;
      }

      if (SphereUtilities::IsPointOnSurfaceAndMovingOut<Real_v, false>(sphere, point, direction)) {
        distance = Real_v(0.0);
        return;
      }
    }

    Real_v sd1(kInfLength);
    Real_v sd2(kInfLength);
    Real_v d2 = (pDotV3d * pDotV3d - c);
    if (d2 < Real_v(0.) || ((c > Real_v(0.)) && (pDotV3d > Real_v(0.)))) return;

    sd1 = -pDotV3d - Sqrt(d2);

    Real_v outerDist(kInfLength);
    Real_v innerDist(kInfLength);

    if (sphere.fFullSphere) {
      if (sd1 > Real_v(kTolerance)) outerDist = sd1;
    } else {
      if (sd1 > Real_v(kTolerance) && sd1 < Real_v(kInfLength)) {
        Vector3D<Real_v> tmpPt = point + sd1 * direction;
        // Only active angular cuts can reject a radial entry candidate.
        if ((fullPhiSphere || sphere.fPhiWedge.Inside<Real_v>(tmpPt) != EInside::kOutside) &&
            (fullThetaSphere || sphere.fThetaCone.Inside<Real_v>(tmpPt) != EInside::kOutside))
          outerDist = sd1;
      }
    }

    if (sphere.fRmin) {
      c  = rad2 - sphere.fRmin * sphere.fRmin;
      d2 = pDotV3d * pDotV3d - c;
      if (d2 >= Real_v(0.)) sd2 = -pDotV3d + Sqrt(d2);
      const bool innerTangentContinuation = (c <= Real_v(2. * kTolerance) * sphere.fRmin) &&
                                            (c >= -Real_v(2. * kTolerance) * sphere.fRmin) &&
                                            (pDotV3d >= -kToleranceDist<Real_v>);

      if (sphere.fFullSphere) {
        // Exact tangents on a concave inner sphere enter material at second
        // order; DistanceToOut owns the complementary hollow-crossing side.
        if (sd2 > Real_v(kTolerance) || (innerTangentContinuation && sd2 >= -Real_v(kTolerance))) {
          innerDist = innerTangentContinuation ? Real_v(0.) : sd2;
        }
      } else {
        if (innerTangentContinuation || (sd2 > Real_v(kTolerance) && sd2 < Real_v(kInfLength))) {
          const Real_v innerCandidate = innerTangentContinuation ? Real_v(0.) : sd2;
          Vector3D<Real_v> tmpPt      = point + innerCandidate * direction;
          if ((fullPhiSphere || sphere.fPhiWedge.Inside<Real_v>(tmpPt) != EInside::kOutside) &&
              (fullThetaSphere || sphere.fThetaCone.Inside<Real_v>(tmpPt) != EInside::kOutside))
            innerDist = innerCandidate;
        }
      }
    }

    distance = Min(outerDist, innerDist);

    if (!fullPhiSphere) {
      GetMinDistFromPhi<Real_v, true>(sphere, point, direction, done, distance);
    }

    Real_v distThetaMin(kInfLength);

    if (!fullThetaSphere) {
      bool intsect1 = false;
      bool intsect2 = false;
      Real_v distTheta1(kInfLength), distCone1(kInfLength);
      Real_v distTheta2(kInfLength), distCone2(kInfLength);
      Real_v rmin2 = sphere.fRmin * sphere.fRmin;
      Real_v rmax2 = sphere.fRmax * sphere.fRmax;
      Vector3D<Real_v> coneIntSecPt1;
      Vector3D<Real_v> coneIntSecPt2;

      sphere.fThetaCone.DistanceToIn<Real_v>(point, direction, distTheta1, distTheta2, intsect1,
                                             intsect2); //,cone1IntSecPt, cone2IntSecPt);
      bool cone1MovesIn    = false;
      bool cone1EntersApex = false;
      bool cone1ZeroRoot   = intsect1 && (Abs(distTheta1) <= Real_v(kTolerance));
      bool candidateCone1  = intsect1 && ((distTheta1 > Real_v(kTolerance)) || cone1ZeroRoot);
      Real_v pDotV2d       = point.x() * direction.x() + point.y() * direction.y();
      Real_v dirRho2       = direction.Perp2();
      Real_v pointZDir     = point.z() * direction.z();
      Real_v dirZ2         = direction.z() * direction.z();

      if (candidateCone1) {
        distCone1 = rad2 + distTheta1 * (Real_v(2.) * pDotV3d + distTheta1);
        candidateCone1 &= distCone1 >= rmin2 && distCone1 <= rmax2;
        if (candidateCone1) {
          Real_v motion1 = -direction.z();
          Real_v a       = Real_v(0.);
          if (Abs(sphere.fSTheta - kHalfPi) > kTolerance) {
            Real_v tanTheta2 = Real_v(sphere.fThetaCone.GetTanSTheta2());
            a                = dirRho2 - dirZ2 * tanTheta2;
            Real_v b         = pDotV2d - pointZDir * tanTheta2;
            motion1          = b + distTheta1 * a;
          }
          cone1MovesIn = SphereUtilities::IsThetaConeMotion<Real_v, true, false>(sphere, motion1);
          // A zero theta-root with unresolved first-order motion can still
          // enter material through the cone curvature.
          if (cone1ZeroRoot && !cone1MovesIn && Abs(motion1) <= kToleranceDist<Real_v>) {
            cone1MovesIn = (sphere.fSTheta <= kHalfPi) ? (a > Real_v(0.)) : (a < Real_v(0.));
          }
          cone1EntersApex = (sphere.fRmin == 0.) && (distCone1 <= Real_v(kTolerance * kTolerance)) &&
                            sphere.fThetaCone.Contains<Real_v>(direction);
          if (!fullPhiSphere) cone1EntersApex &= sphere.fPhiWedge.Contains<Real_v>(direction);
        }
      }

      bool cone2MovesIn    = false;
      bool cone2EntersApex = false;
      bool cone2ZeroRoot   = intsect2 && (Abs(distTheta2) <= Real_v(kTolerance));
      bool candidateCone2  = intsect2 && ((distTheta2 > Real_v(kTolerance)) || cone2ZeroRoot);
      if (candidateCone2) {
        distCone2 = rad2 + distTheta2 * (Real_v(2.) * pDotV3d + distTheta2);
        candidateCone2 &= distCone2 >= rmin2 && distCone2 <= rmax2;
        if (candidateCone2) {
          Real_v motion2 = -direction.z();
          Real_v a       = Real_v(0.);
          if (Abs(sphere.eTheta - kHalfPi) > kTolerance) {
            Real_v tanTheta2 = Real_v(sphere.fThetaCone.GetTanETheta2());
            a                = dirRho2 - dirZ2 * tanTheta2;
            Real_v b         = pDotV2d - pointZDir * tanTheta2;
            motion2          = b + distTheta2 * a;
          }
          cone2MovesIn = SphereUtilities::IsThetaConeMotion<Real_v, false, false>(sphere, motion2);
          // Same zero-root rule for the end-theta cone.
          if (cone2ZeroRoot && !cone2MovesIn && Abs(motion2) <= kToleranceDist<Real_v>) {
            cone2MovesIn = (sphere.eTheta <= kHalfPi) ? (a < Real_v(0.)) : (a > Real_v(0.));
          }
          cone2EntersApex = (sphere.fRmin == 0.) && (distCone2 <= Real_v(kTolerance * kTolerance)) &&
                            sphere.fThetaCone.Contains<Real_v>(direction);
          if (!fullPhiSphere) cone2EntersApex &= sphere.fPhiWedge.Contains<Real_v>(direction);
        }
      }

      bool isValidCone1 = candidateCone1 && (cone1MovesIn || cone1EntersApex);
      bool isValidCone2 = candidateCone2 && (cone2MovesIn || cone2EntersApex);

      if (!fullPhiSphere && isValidCone1 && !cone1EntersApex) {
        coneIntSecPt1 = point + distTheta1 * direction;
        isValidCone1  = sphere.fPhiWedge.Contains<Real_v>(coneIntSecPt1);
      }
      if (!fullPhiSphere && isValidCone2 && !cone2EntersApex) {
        coneIntSecPt2 = point + distTheta2 * direction;
        isValidCone2  = sphere.fPhiWedge.Contains<Real_v>(coneIntSecPt2);
      }
      if (isValidCone2 && !isValidCone1) distThetaMin = distTheta2;
      if (isValidCone1 && !isValidCone2) distThetaMin = distTheta1;
      if (isValidCone1 && isValidCone2) distThetaMin = Min(distTheta1, distTheta2);

      bool isApexEntry = (sphere.fRmin == 0.) && (rad2 <= Real_v(kTolerance * kTolerance)) &&
                         sphere.fThetaCone.Contains<Real_v>(direction);
      if (!fullPhiSphere) isApexEntry &= sphere.fPhiWedge.Contains<Real_v>(direction);
      if (isApexEntry) distThetaMin = Real_v(0.);
    }

    distance = Min(distThetaMin, distance);
  }

  /// @brief Compute the first valid exit distance.
  /// @details Considers radial shell roots, theta-cone roots, and phi-plane
  /// roots, returning the nearest accepted boundary crossing.
  /// @tparam RecordHitInfo Whether to populate the optional exit information.
  /// @param sphere Sphere parameters and cached angular state.
  /// @param point Local starting point.
  /// @param direction Normalized local direction.
  /// @param[out] distance Set to the exit distance, `-1`, `0`, or `kInfLength`.
  /// @param hitInfo Optional-output view; required when `RecordHitInfo` is true.
  template <typename Real_v, bool RecordHitInfo>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void DistanceToOutImpl(
      UnplacedStruct_t const &sphere, Vector3D<Real_v> const &point, Vector3D<Real_v> const &direction,
      Real_v const & /* stepMax */, Real_v &distance, SurfaceHitView<Real_v> *hitInfo)
  {

    distance                                = kInfLength;
    [[maybe_unused]] ExitSurface hitSurface = ExitSurface::kNone;
    if constexpr (RecordHitInfo) hitInfo->Clear();

    Real_v snxt(kInfLength);

    Real_v d2(0.);

    Real_v pDotV3d = point.Dot(direction);

    Real_v rad2 = point.Mag2();
    Real_v c    = rad2 - sphere.fRmax * sphere.fRmax;

    Real_v sd1(kInfLength);
    Real_v sd2(kInfLength);

    const Real_v outerLimit = sphere.fRmax + Real_v(kTolerance);
    bool completelyOutside  = rad2 > outerLimit * outerLimit;
    if (sphere.fRmin > kTolerance) {
      const Real_v innerLimit = sphere.fRmin - Real_v(kTolerance);
      completelyOutside |= rad2 < innerLimit * innerLimit;
    }
    if (!sphere.fFullPhiSphere) {
      bool unused = false, outsidePhi = false;
      sphere.fPhiWedge.GenericKernelForContainsAndInside<Real_v, false>(point, unused, outsidePhi);
      if (outsidePhi) {
        bool insidePhi = false;
        sphere.fPhiWedge.GenericKernelForContainsAndInside<Real_v, true>(point, insidePhi, outsidePhi);
      }
      completelyOutside |= outsidePhi;
    }
    bool thetaNeedsSurfaceCheck = false;
    if (!sphere.fFullThetaSphere) {
      const Real_v radius = Sqrt(rad2);
      // Match ThetaCone's 2*10*kAngTolerance broad band before invoking its exact surface classifier.
      const Real_v margin = Real_v(20. * sphere.kAngTolerance);
      const bool clearlyInsideTheta =
          point.z() < radius * sphere.cosSTheta - margin && point.z() > radius * sphere.cosETheta + margin;
      thetaNeedsSurfaceCheck = !clearlyInsideTheta;
      if (thetaNeedsSurfaceCheck) completelyOutside |= sphere.fThetaCone.IsCompletelyOutside<Real_v>(point);
    }
    if (completelyOutside) {
      distance = Real_v(-1.0);
      return;
    }

    // Preserve zero exits throughout the radial surface band; angular containment is needed only there.
    const bool outerSurfaceExit =
        rad2 >= MakeMinusTolerantSquare<true>(sphere.fRmax) && pDotV3d > Real_v(-kToleranceDist<Real_v>);
    const bool innerSurfaceExit =
        sphere.fRmin && rad2 <= MakePlusTolerantSquare<true>(sphere.fRmin) && pDotV3d < Real_v(-kToleranceDist<Real_v>);
    if ((outerSurfaceExit || innerSurfaceExit) && (sphere.fFullPhiSphere || sphere.fPhiWedge.Contains<Real_v>(point)) &&
        (sphere.fFullThetaSphere || sphere.fThetaCone.Contains<Real_v>(point))) {
      distance = Real_v(0.);
      if constexpr (RecordHitInfo)
        SetExitHit(sphere, point, direction, distance, outerSurfaceExit ? ExitSurface::kOuter : ExitSurface::kInner,
                   hitInfo);
      return;
    }

    d2 = (pDotV3d * pDotV3d - c);
    if (d2 >= Real_v(0.)) sd1 = -pDotV3d + Sqrt(d2);

    if (sphere.fRmin) {
      c  = rad2 - sphere.fRmin * sphere.fRmin;
      d2 = (pDotV3d * pDotV3d - c);
      if (d2 >= Real_v(0.) && pDotV3d < -kToleranceDist<Real_v>) sd2 = -pDotV3d - Sqrt(d2);
    }

    if constexpr (RecordHitInfo) {
      if (sd1 < sd2) {
        snxt       = sd1;
        hitSurface = ExitSurface::kOuter;
      } else {
        snxt       = sd2;
        hitSurface = sd2 < Real_v(kInfLength) ? ExitSurface::kInner : ExitSurface::kNone;
      }
    } else {
      snxt = Min(sd1, sd2);
    }

    const bool upperHemisphere = sphere.fSTheta == 0. && sphere.eTheta == kHalfPi;
    const bool lowerHemisphere = sphere.fSTheta == kHalfPi && sphere.eTheta == kPi;
    if (sphere.fFullPhiSphere && (upperHemisphere || lowerHemisphere)) {
      distance              = snxt;
      const bool exitsPlane = upperHemisphere ? direction.z() < Real_v(0.) : direction.z() > Real_v(0.);
      if (exitsPlane) {
        Real_v planeDistance = -point.z() / direction.z();
        if (Abs(point.z()) < Real_v(kTolerance)) planeDistance = Real_v(0.);
        if (planeDistance < distance) {
          distance = planeDistance;
          if constexpr (RecordHitInfo) hitSurface = upperHemisphere ? ExitSurface::kThetaEnd : ExitSurface::kThetaStart;
        }
      }
      if constexpr (RecordHitInfo) SetExitHit(sphere, point, direction, distance, hitSurface, hitInfo);
      return;
    }

    Real_v distThetaMin(kInfLength);
    Real_v distPhiMin(kInfLength);
    [[maybe_unused]] ExitSurface thetaSurface = ExitSurface::kNone;

    if (!sphere.fFullThetaSphere) {
      bool intsect1 = false;
      bool intsect2 = false;
      Real_v distTheta1(kInfLength);
      Real_v distTheta2(kInfLength);
      ThetaDistancesToOut(sphere, point, direction, distTheta1, distTheta2, intsect1, intsect2);
      if (thetaNeedsSurfaceCheck) {
        if (SphereUtilities::IsPointOnThetaSurfaceAndMovingOut<Real_v, true, true>(sphere, point, direction)) {
          distTheta1 = Real_v(0.);
          intsect1   = true;
        }
        if (SphereUtilities::IsPointOnThetaSurfaceAndMovingOut<Real_v, false, true>(sphere, point, direction)) {
          distTheta2 = Real_v(0.);
          intsect2   = true;
        }
      }
      // The cone helper can expose a backward root; only tolerated surface roots belong to this exit query.
      intsect1 &= distTheta1 >= -Real_v(kTolerance);
      intsect2 &= distTheta2 >= -Real_v(kTolerance);
      if (intsect1 && Abs(distTheta1) <= Real_v(kTolerance)) {
        const bool apexExit = rad2 <= Real_v(kTolerance * kTolerance) && direction.z() > Real_v(sphere.cosSTheta);
        Real_v motion1      = -direction.z();
        Real_v a            = Real_v(0.);
        if (Abs(sphere.fSTheta - kHalfPi) > kTolerance) {
          Real_v tanTheta2 = Real_v(sphere.fThetaCone.GetTanSTheta2());
          a                = direction.Perp2() - direction.z() * direction.z() * tanTheta2;
          motion1 = point.x() * direction.x() + point.y() * direction.y() - point.z() * direction.z() * tanTheta2;
        }
        // From a theta surface, an inward/tangential material continuation
        // must exit at the next boundary, not at the current zero root.
        bool materialContinuation = SphereUtilities::IsThetaConeMotion<Real_v, true, false>(sphere, motion1);
        if (!materialContinuation && Abs(motion1) <= kToleranceDist<Real_v>) {
          materialContinuation = (sphere.fSTheta <= kHalfPi) ? (a > Real_v(0.)) : (a < Real_v(0.));
        }
        if (materialContinuation && !apexExit) intsect1 = false;
      }
      if (intsect2 && Abs(distTheta2) <= Real_v(kTolerance)) {
        const bool apexExit = rad2 <= Real_v(kTolerance * kTolerance) && direction.z() < Real_v(sphere.cosETheta);
        Real_v motion2      = -direction.z();
        Real_v a            = Real_v(0.);
        if (Abs(sphere.eTheta - kHalfPi) > kTolerance) {
          Real_v tanTheta2 = Real_v(sphere.fThetaCone.GetTanETheta2());
          a                = direction.Perp2() - direction.z() * direction.z() * tanTheta2;
          motion2 = point.x() * direction.x() + point.y() * direction.y() - point.z() * direction.z() * tanTheta2;
        }
        // Same zero-root suppression for end-theta material continuation.
        bool materialContinuation = SphereUtilities::IsThetaConeMotion<Real_v, false, false>(sphere, motion2);
        if (!materialContinuation && Abs(motion2) <= kToleranceDist<Real_v>) {
          materialContinuation = (sphere.eTheta <= kHalfPi) ? (a < Real_v(0.)) : (a > Real_v(0.));
        }
        if (materialContinuation && !apexExit) intsect2 = false;
      }
      if (intsect2 && !intsect1) {
        distThetaMin = distTheta2;
        if constexpr (RecordHitInfo) thetaSurface = ExitSurface::kThetaEnd;
      }
      if (!intsect2 && intsect1) {
        distThetaMin = distTheta1;
        if constexpr (RecordHitInfo) thetaSurface = ExitSurface::kThetaStart;
      }
      if (intsect2 && intsect1) {
        distThetaMin = Min(distTheta1, distTheta2);
        if constexpr (RecordHitInfo)
          thetaSurface = distTheta1 < distTheta2 ? ExitSurface::kThetaStart : ExitSurface::kThetaEnd;
      }
    }

    if constexpr (RecordHitInfo) {
      distance = snxt;
      if (distThetaMin < distance) {
        distance   = distThetaMin;
        hitSurface = thetaSurface;
      }
    } else {
      distance = Min(distThetaMin, snxt);
    }

    if (!sphere.fFullPhiSphere) {
      if (sphere.fDPhi <= kPi) {
        Real_v distPhi1;
        Real_v distPhi2;
        PhiDistancesToOut(sphere, point, direction, distPhi1, distPhi2);
        distPhiMin = Min(distPhi1, distPhi2);
        if constexpr (RecordHitInfo) {
          if (distPhiMin < distance) {
            distance   = distPhiMin;
            hitSurface = distPhi1 < distPhi2 ? ExitSurface::kPhiStart : ExitSurface::kPhiEnd;
          }
        } else {
          distance = Min(distPhiMin, distance);
        }
      } else {
        GetMinDistToOutFromPhi<Real_v, RecordHitInfo>(sphere, point, direction, distance, &hitSurface);
      }
    }
    // A tolerance-surface ray can have no quadratic root when its direction is nearly tangent. Keep the complete
    // surface classification as a cold fallback rather than paying for it before every regular distance query.
    if (distance >= Real_v(kInfLength) &&
        SphereUtilities::IsPointOnSurfaceAndMovingOut<Real_v, true>(sphere, point, direction)) {
      distance = Real_v(0.);
      if constexpr (RecordHitInfo) hitSurface = SelectSurfaceExit(sphere, point, direction);
    }
    if constexpr (RecordHitInfo) SetExitHit(sphere, point, direction, distance, hitSurface, hitInfo);
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void DistanceToOutWithHit(UnplacedStruct_t const &sphere,
                                                                                Vector3D<Real_v> const &point,
                                                                                Vector3D<Real_v> const &direction,
                                                                                Real_v const &stepMax, Real_v &distance,
                                                                                SurfaceHitView<Real_v> *hitInfo)
  {
    DistanceToOutImpl<Real_v, true>(sphere, point, direction, stepMax, distance, hitInfo);
  }

  /// @brief Compute the first valid exit distance and optionally return a supporting normal.
  /// @param sphere Sphere parameters and cached angular state.
  /// @param point Local starting point.
  /// @param direction Normalized local direction.
  /// @param stepMax Currently unused distance limit.
  /// @param distance Set to the exit distance, `-1`, `0`, or `kInfLength`.
  /// @param hitInfo Optional view receiving `leftBehind` and, when requested and
  /// valid, the outward unit normal. Passing null preserves the ordinary query.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void DistanceToOut(UnplacedStruct_t const &sphere,
                                                                         Vector3D<Real_v> const &point,
                                                                         Vector3D<Real_v> const &direction,
                                                                         Real_v const &stepMax, Real_v &distance,
                                                                         SurfaceHitView<Real_v> *hitInfo = nullptr)
  {
    if (hitInfo)
      DistanceToOutWithHit(sphere, point, direction, stepMax, distance, hitInfo);
    else
      DistanceToOutImpl<Real_v, false>(sphere, point, direction, stepMax, distance, nullptr);
  }

  /// Test the two half-planes of a wide phi cut, stopping once a candidate cannot beat another exit.
  /// The current distance is already the nearest radial/theta exit. A phi candidate preceding it therefore remains
  /// inside those extents and needs only the selected half-plane check.
  template <typename Real_v, bool RecordHitInfo>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void GetMinDistToOutFromPhi(UnplacedStruct_t const &sphere,
                                                                                  Vector3D<Real_v> const &localPoint,
                                                                                  Vector3D<Real_v> const &localDir,
                                                                                  Real_v &distance,
                                                                                  ExitSurface *hitSurface)
  {
    Real_v distPhi1(kInfLength);
    Real_v distPhi2(kInfLength);
    PhiDistancesToOut(sphere, localPoint, localDir, distPhi1, distPhi2);

    const Real_v nearest  = Min(distPhi1, distPhi2);
    const Real_v farthest = Max(distPhi1, distPhi2);
    if (nearest < distance) {
      const Vector3D<Real_v> hit = localPoint + nearest * localDir;
      const bool onFace          = (nearest == distPhi1) ? sphere.fPhiWedge.IsOnSurfaceGeneric<Real_v, true>(hit)
                                                         : sphere.fPhiWedge.IsOnSurfaceGeneric<Real_v, false>(hit);
      if (onFace) {
        // At the shared axis, leaving one half-space need not leave the wide wedge.
        // Motion into the other half-space excludes this and any later phi exit.
        if (hit.Perp2() < Real_v(kTolerance * kTolerance)) {
          const auto other = nearest == distPhi1 ? sphere.fPhiWedge.GetNormal2() : sphere.fPhiWedge.GetNormal1();
          if (localDir.x() * other.x() + localDir.y() * other.y() > Real_v(1.e-12)) return;
        }
        distance = nearest;
        if constexpr (RecordHitInfo) *hitSurface = distPhi1 < distPhi2 ? ExitSurface::kPhiStart : ExitSurface::kPhiEnd;
        return;
      }
    }
    if (farthest < distance) {
      const Vector3D<Real_v> hit = localPoint + farthest * localDir;
      const bool onFace          = (farthest == distPhi1) ? sphere.fPhiWedge.IsOnSurfaceGeneric<Real_v, true>(hit)
                                                          : sphere.fPhiWedge.IsOnSurfaceGeneric<Real_v, false>(hit);
      if (onFace) {
        distance = farthest;
        if constexpr (RecordHitInfo) *hitSurface = distPhi1 > distPhi2 ? ExitSurface::kPhiStart : ExitSurface::kPhiEnd;
      }
    }
  }

  /// @brief Test phi-plane crossings and update the best distance.
  /// @details Candidate phi intersections are accepted only when the hit point
  /// is on the selected phi surface and inside the radial/theta extents.
  /// @tparam DistToIn Selects entry or exit phi-plane queries.
  /// @param sphere Sphere parameters and cached angular state.
  /// @param localPoint Local starting point.
  /// @param localDir Normalized local direction.
  /// @param done Suppresses updates when an earlier phase completed the query.
  /// @param distance Current best distance, updated if a closer phi hit exists.
  /// @param hitSurface Winning phi boundary, updated only when hit recording is enabled.
  template <typename Real_v, bool DistToIn, bool RecordHitInfo = false>
  VECCORE_ATT_HOST_DEVICE static void GetMinDistFromPhi(UnplacedStruct_t const &sphere,
                                                        Vector3D<Real_v> const &localPoint,
                                                        Vector3D<Real_v> const &localDir, bool done, Real_v &distance,
                                                        ExitSurface *hitSurface = nullptr)
  {
    Real_v distPhi1(kInfLength);
    Real_v distPhi2(kInfLength);
    Real_v dist(kInfLength);

    if (DistToIn)
      sphere.fPhiWedge.DistanceToIn<Real_v>(localPoint, localDir, distPhi1, distPhi2);
    else
      sphere.fPhiWedge.DistanceToOut<Real_v>(localPoint, localDir, distPhi1, distPhi2);

    bool containsCond1 = false, containsCond2 = false;
    bool tempCond = false;
    // Min Face
    dist = Min(distPhi1, distPhi2);
    if (dist < kInfLength) {
      Vector3D<Real_v> tmpPt = localPoint + dist * localDir;
      Real_v rad2            = tmpPt.Mag2();

      tempCond = ((dist == distPhi1) && sphere.fPhiWedge.IsOnSurfaceGeneric<Real_v, true>(tmpPt)) ||
                 ((dist == distPhi2) && sphere.fPhiWedge.IsOnSurfaceGeneric<Real_v, false>(tmpPt));

      containsCond1 = tempCond && (rad2 > sphere.fRmin * sphere.fRmin) && (rad2 < sphere.fRmax * sphere.fRmax) &&
                      sphere.fThetaCone.Contains<Real_v>(tmpPt);
    }

    if (!done && containsCond1) {
      if constexpr (RecordHitInfo) {
        if (dist < distance) {
          distance    = dist;
          *hitSurface = distPhi1 < distPhi2 ? ExitSurface::kPhiStart : ExitSurface::kPhiEnd;
        }
      } else {
        distance = Min(dist, distance);
      }
    }

    // Max Face
    dist = Max(distPhi1, distPhi2);
    if (dist < kInfLength) {
      Vector3D<Real_v> tmpPt = localPoint + dist * localDir;

      Real_v rad2 = tmpPt.Mag2();
      tempCond    = ((dist == distPhi1) && sphere.fPhiWedge.IsOnSurfaceGeneric<Real_v, true>(tmpPt)) ||
                    ((dist == distPhi2) && sphere.fPhiWedge.IsOnSurfaceGeneric<Real_v, false>(tmpPt));

      containsCond2 = tempCond && (rad2 > sphere.fRmin * sphere.fRmin) && (rad2 < sphere.fRmax * sphere.fRmax) &&
                      sphere.fThetaCone.Contains<Real_v>(tmpPt);
    }
    if (!done && !containsCond1 && containsCond2) {
      if constexpr (RecordHitInfo) {
        if (dist < distance) {
          distance    = dist;
          *hitSurface = distPhi1 > distPhi2 ? ExitSurface::kPhiStart : ExitSurface::kPhiEnd;
        }
      } else {
        distance = Min(dist, distance);
      }
    }
  }

  /// @brief Compute signed safety for entering the sphere.
  /// @details Full-theta spheres use radial safety, optionally combined with
  /// phi safety, to decide outside, surface, and inside states directly.
  /// @param sphere Sphere parameters and cached angular state.
  /// @param point Local point to test.
  /// @param safety Set to positive distance, `0`, or `-1`.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SafetyToIn(UnplacedStruct_t const &sphere,
                                                                      Vector3D<Real_v> const &point, Real_v &safety)
  {
    // General Precalcs
    Real_v rad = point.Mag();

    Real_v safeRMin(0.);
    Real_v safeRMax(0.);

    if (sphere.fFullThetaSphere) {
      // With no theta cuts, signed radial/phi safety fully determines inside,
      // surface, and outside state, so avoid a separate classification pass.
      if (sphere.fRmin) {
        safeRMin = sphere.fRmin - rad;
        safeRMax = rad - sphere.fRmax;
        safety   = safeRMin > safeRMax ? safeRMin : safeRMax;
      } else {
        safety = rad - sphere.fRmax;
      }

      if (!sphere.fFullPhiSphere) safety = Max(sphere.fPhiWedge.SafetyToIn<Real_v>(point), safety);
      if (safety > Real_v(kTolerance)) return;
      safety = safety < Real_v(-kTolerance) ? Real_v(-1.) : Real_v(0.);
      return;
    }

    bool completelyinside  = false;
    bool completelyoutside = false;
    GenericKernelForContainsAndInside<Real_v, true>(sphere, point, completelyinside, completelyoutside);

    if (completelyinside) {
      safety = Real_v(-1.0);
      return;
    }

    if (!completelyoutside) {
      safety = Real_v(0.0);
      return;
    }

    if (sphere.fRmin) {
      safeRMin = sphere.fRmin - rad;
      safeRMax = rad - sphere.fRmax;
      safety   = safeRMin > safeRMax ? safeRMin : safeRMax;
    } else {
      safety = rad - sphere.fRmax;
    }
    // Distance to r shells over

    // Distance to phi extent
    if (!sphere.fFullPhiSphere) {
      Real_v safetyPhi = sphere.fPhiWedge.SafetyToIn<Real_v>(point);
      safety           = Max(safetyPhi, safety);
    }

    // Distance to Theta extent
    if (!sphere.fFullThetaSphere) {
      Real_v safetyTheta = sphere.fThetaCone.SafetyToIn<Real_v>(point);
      safety             = Max(safetyTheta, safety);
    }
  }

  /// @brief Compute signed safety for leaving the sphere.
  /// @details Uses the nearest active radial, phi, or theta boundary for points
  /// classified safely inside.
  /// @param sphere Sphere parameters and cached angular state.
  /// @param point Local point to test.
  /// @param safety Set to positive distance, `0`, or `-1`.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SafetyToOut(UnplacedStruct_t const &sphere,
                                                                       Vector3D<Real_v> const &point, Real_v &safety)
  {

    Real_v rad = point.Mag();

    bool completelyinside  = false;
    bool completelyoutside = false;
    GenericKernelForContainsAndInside<Real_v, true>(sphere, point, completelyinside, completelyoutside);
    if (completelyoutside) {
      safety = Real_v(-1.0);
      return;
    }

    if (!completelyinside) {
      safety = Real_v(0.0);
      return;
    }

    // Distance to r shells
    if (sphere.fRmin) {
      Real_v safeRMin = (rad - sphere.fRmin);
      Real_v safeRMax = (sphere.fRmax - rad);
      safety          = safeRMin < safeRMax ? safeRMin : safeRMax;
    } else {
      safety = sphere.fRmax - rad;
    }

    // Distance to phi extent
    if (!sphere.fFullPhiSphere) {
      Real_v safetyPhi = sphere.fPhiWedge.SafetyToOut<Real_v>(point);
      safety           = Min(safetyPhi, safety);
    }

    // Distance to Theta extent
    if (!sphere.fFullThetaSphere) {
      Real_v safeTheta = sphere.fThetaCone.SafetyToOut<Real_v>(point);
      safety           = Min(safeTheta, safety);
    }
  }

  /// @brief Approximate the normal from the closest sphere boundary.
  /// @details Used by `Normal` when the point is not on a recognized surface.
  /// @param sphere Sphere parameters and cached angular state.
  /// @param point Local point to evaluate.
  /// @return Normal of the closest radial, phi, or theta boundary.
  template <typename Real_v>
  VECCORE_ATT_HOST_DEVICE static Vector3D<Real_v> ApproxSurfaceNormalKernel(UnplacedStruct_t const &sphere,
                                                                            Vector3D<Real_v> const &point)
  {
    Vector3D<Real_v> norm(0., 0., 0.);
    Real_v radius     = point.Mag();
    Real_v distRMax   = Abs(radius - sphere.fRmax);
    Real_v distRMin   = InfinityLength<Real_v>();
    Real_v distPhi1   = InfinityLength<Real_v>();
    Real_v distPhi2   = InfinityLength<Real_v>();
    Real_v distTheta1 = InfinityLength<Real_v>();
    Real_v distTheta2 = InfinityLength<Real_v>();
    Real_v distMin    = distRMax;

    if (sphere.fRmin > 0.) {
      distRMin = Abs(sphere.fRmin - radius);
      distMin  = Min(distRMin, distRMax);
    }

    if (!sphere.fFullPhiSphere) {
      distPhi1 = Abs(point.x() * sphere.fPhiWedge.GetNormal1().x() + point.y() * sphere.fPhiWedge.GetNormal1().y());
      distPhi2 = Abs(point.x() * sphere.fPhiWedge.GetNormal2().x() + point.y() * sphere.fPhiWedge.GetNormal2().y());
      distMin  = Min(distMin, distPhi1, distPhi2);
    }

    if (!sphere.fFullThetaSphere) {
      Real_v rho = point.Perp();
      distTheta1 = sphere.fThetaCone.DistanceToLine<Real_v>(sphere.fThetaCone.GetSlope1(), rho, point.z());
      distTheta2 = sphere.fThetaCone.DistanceToLine<Real_v>(sphere.fThetaCone.GetSlope2(), rho, point.z());
      distMin    = Min(distMin, distTheta1, distTheta2);
    }

    if (distMin == distRMax) norm = point.Unit();
    if (distMin == distRMin) norm = -point.Unit();

    Vector3D<Real_v> normal1 = sphere.fPhiWedge.GetNormal1();
    Vector3D<Real_v> normal2 = sphere.fPhiWedge.GetNormal2();
    if (distMin == distPhi1) norm = -normal1;
    if (distMin == distPhi2) norm = -normal2;

    if (distMin == distTheta1) norm += sphere.fThetaCone.GetNormal1<Real_v>(point);
    if (distMin == distTheta2) norm += sphere.fThetaCone.GetNormal2<Real_v>(point);

    return norm;
  }

  /// @brief Return the outward normal and whether a surface was identified.
  /// @details Surface normals are accumulated from all matching active
  /// boundaries; off-surface points use the closest-boundary fallback.
  /// @param sphere Sphere parameters and cached angular state.
  /// @param point Local surface point.
  /// @param valid Set when at least one surface was identified.
  /// @return Outward normal direction.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static Vector3D<Real_v> Normal(UnplacedStruct_t const &sphere,
                                                                              Vector3D<Real_v> const &point,
                                                                              bool &valid)
  {
    Vector3D<Real_v> normal(0., 0., 0.);
    normal.Set(1e-30);

    /* Assumption : This function assumes that the point is on the surface.
     *
     * Algorithm :
     * Detect all those surfaces on which the point is at, and count the
     * numOfSurfaces. if(numOfSurfaces == 1) then normal corresponds to the
     * normal for that particular case.
     *
     * if(numOfSurfaces > 1 ), then add the normals corresponds to different
     * cases, and finally normalize it and return.
     *
     * We need following function
     * IsPointOnInnerRadius()
     * IsPointOnOuterRadius()
     * IsPointOnStartPhi()
     * IsPointOnEndPhi()
     * IsPointOnStartTheta()
     * IsPointOnEndTheta()
     *
     * set valid=true if numOfSurface > 0
     *
     * if above mentioned assumption not followed , ie.
     * In case the given point is outside, then find the closest boundary,
     * the required normal will be the normal to that boundary.
     * This logic is implemented in "ApproxSurfaceNormalKernel" function
     */

    bool isPointOutside = false;

    bool isPointInside = false;
    GenericKernelForContainsAndInside<Real_v, true>(sphere, point, isPointInside, isPointOutside);
    if (isPointInside || isPointOutside) normal = ApproxSurfaceNormalKernel<Real_v>(sphere, point);
    bool needsNormalization = isPointInside || isPointOutside;

    valid = false;

    int noSurfaces            = 0;
    bool isPointOnOuterRadius = SphereUtilities::IsPointOnOuterRadius<Real_v>(sphere, point);

    if (isPointOnOuterRadius) ++noSurfaces;
    if (!isPointOutside && isPointOnOuterRadius) normal += point.Unit();

    if (sphere.fRmin) {
      bool isPointOnInnerRadius = SphereUtilities::IsPointOnInnerRadius<Real_v>(sphere, point);
      if (isPointOnInnerRadius) ++noSurfaces;
      if (!isPointOutside && isPointOnInnerRadius) normal -= point.Unit();
    }

    if (!sphere.fFullPhiSphere) {
      bool isPointOnStartPhi = SphereUtilities::IsPointOnStartPhi<Real_v>(sphere, point);
      bool isPointOnEndPhi   = SphereUtilities::IsPointOnEndPhi<Real_v>(sphere, point);
      if (isPointOnStartPhi) ++noSurfaces;
      if (isPointOnEndPhi) ++noSurfaces;
      if (!isPointOutside && isPointOnStartPhi) normal -= sphere.fPhiWedge.GetNormal1();
      if (!isPointOutside && isPointOnEndPhi) normal -= sphere.fPhiWedge.GetNormal2();
    }

    if (!sphere.fFullThetaSphere) {
      bool isPointOnStartTheta = SphereUtilities::IsPointOnStartTheta<Real_v>(sphere, point);
      bool isPointOnEndTheta   = SphereUtilities::IsPointOnEndTheta<Real_v>(sphere, point);
      needsNormalization |= isPointOnStartTheta || isPointOnEndTheta;

      if (isPointOnStartTheta) ++noSurfaces;
      if (!isPointOutside && isPointOnStartTheta) normal += sphere.fThetaCone.GetNormal1<Real_v>(point);

      if (isPointOnEndTheta) ++noSurfaces;
      if (!isPointOutside && isPointOnEndTheta) normal += sphere.fThetaCone.GetNormal2<Real_v>(point);

      Vector3D<Real_v> tempNormal(0., 0., -1.);
      if (!isPointOutside && isPointOnStartTheta && isPointOnEndTheta && (sphere.eTheta <= kPi / 2.))
        normal = tempNormal;
      Vector3D<Real_v> tempNormal2(0., 0., 1.);
      if (!isPointOutside && isPointOnStartTheta && isPointOnEndTheta && (sphere.fSTheta >= kPi / 2.))
        normal = tempNormal2;
    }

    // A unique radial or phi contribution is already unit length.
    if (needsNormalization || noSurfaces != 1) normal.Normalize();

    valid = (noSurfaces > 0);

    return normal;
  }
};
} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif // VECGEOM_VOLUMES_KERNEL_sphereIMPLEMENTATION_H_
