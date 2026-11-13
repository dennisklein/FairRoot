/*
 * BooleanImplementation.h
 */

#ifndef BOOLEANINTERSECTIONIMPLEMENTATION_H_
#define BOOLEANINTERSECTIONIMPLEMENTATION_H_

#include "VecGeom/base/Global.h"
#include "VecGeom/base/Vector3D.h"
#include "VecGeom/volumes/BooleanStruct.h"

namespace vecgeom {

inline namespace VECGEOM_IMPL_NAMESPACE {

/**
 * @brief Kernel implementation for Boolean intersection volumes.
 * @details The implementation delegates to the two placed constituents and
 * combines their answers using the intersection convention: a point belongs to
 * the solid only when it belongs to both constituents.
 */
template <>
struct BooleanImplementation<kIntersection> {
  using PlacedShape_t    = PlacedBooleanVolume<kIntersection>;
  using UnplacedVolume_t = UnplacedBooleanVolume<kIntersection>;
  using UnplacedStruct_t = BooleanStruct;

  /**
   * @brief Test whether a point is contained in both constituents.
   * @details Contains is boundary-inclusive for the constituent queries, so a
   * point on a constituent boundary can still be contained by the intersection.
   */
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void Contains(BooleanStruct const &unplaced,
                                                                    Vector3D<Real_v> const &point, bool &inside)
  {
    const bool insideA = unplaced.fLeftVolume->Contains(point);
    const bool insideB = unplaced.fRightVolume->Contains(point);
    inside             = insideA && insideB;
  }

  /**
   * @brief Classify a point with the Boolean intersection inside convention.
   * @details A point outside either constituent is outside the intersection.
   * It is strictly inside only if both constituents classify it as inside; any
   * inside/surface combination that is not outside is on the intersection
   * surface. Cached operand scores gate conservative bounds screening before
   * classification.
   */
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void Inside(BooleanStruct const &unplaced,
                                                                  Vector3D<Real_v> const &point,
                                                                  vecgeom::Inside_t &inside)
  {
    // now use the Inside functionality of left and right components
    // algorithm taken from Geant4 implementation
    VPlacedVolume const *const fPtrSolidA = unplaced.fLeftVolume;
    VPlacedVolume const *const fPtrSolidB = unplaced.fRightVolume;

    const auto positionA = BooleanHelper::InsideOperand(fPtrSolidA, unplaced.fLeftScore, point);

    if (positionA == EInside::kOutside) {
      inside = EInside::kOutside;
      return;
    }

    const auto positionB = BooleanHelper::InsideOperand(fPtrSolidB, unplaced.fRightScore, point);
    // The enum order is inside < surface < outside: the larger state
    // gives the intersection classification without enumerating each pair.
    inside = Max(positionA, positionB);
  }

  /**
   * @brief Compute the distance from outside to enter the intersection.
   * @details The boundary-walking algorithm follows alternating constituent
   * entries until the ray is inside both constituents. A strict-inside guard is
   * kept before the walk because DistanceToIn is a wrong-side query for points
   * already inside the intersection, while surface cases are left to the legacy
   * boundary precheck so entering surface rays still return zero. Cached operand
   * scores gate bounds screening for the initial classifications.
   */
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void DistanceToIn(BooleanStruct const &unplaced,
                                                                        Vector3D<Real_v> const &point,
                                                                        Vector3D<Real_v> const &dir,
                                                                        Real_v const &stepMax, Real_v &distance)
  {
    Vector3D<Real_v> hitpoint = point;

    Real_v d1    = 0.;
    Real_v d2    = 0.;
    Real_v snext = 0.0;

    const auto insideLeftState  = BooleanHelper::InsideOperand(unplaced.fLeftVolume, unplaced.fLeftScore, point);
    const auto insideRightState = BooleanHelper::InsideOperand(unplaced.fRightVolume, unplaced.fRightScore, point);
    // For DistanceToIn, constituent surface starts must be re-entered through
    // the normal walking path; otherwise outside approach points on one
    // constituent and inside the other return a false zero for the intersection.
    auto inleft  = insideLeftState == kInside;
    auto inright = insideRightState == kInside;
    if (insideLeftState == kInside && insideRightState == kInside) {
      // Strictly inside both constituents means strictly inside the
      // intersection, so DistanceToIn is called from the wrong side.
      distance = Real_v(-1.0);
      return;
    }

    // just a pre-check before entering main algorithm
    if (inleft && inright) {
      d1 = unplaced.fLeftVolume->PlacedDistanceToOut(hitpoint, dir, stepMax);
      d2 = unplaced.fRightVolume->PlacedDistanceToOut(hitpoint, dir, stepMax);

      // if we are close to a boundary continue
      if (d1 < 2 * kTolerance) inleft = false;
      if (d2 < 2 * kTolerance) inright = false;

      // otherwise exit
      if (inleft && inright) {
        distance = 0.0;
        return;
      }
    }

    // main loop
    while (1) {
      d1 = d2          = 0;
      Real_v boundary1 = 0., boundary2 = 0.;
      if (!inleft) {
        boundary1 = unplaced.fLeftVolume->DistanceToIn(hitpoint, dir);
        // The half-tolerance push ensures progress, but is not part of the
        // physical entry distance. At a zero entry it can round a surface point
        // strictly inside, so retain the unpushed boundary separately.
        d1 = Max(boundary1, Real_v(kHalfTolerance));
        if (d1 > 1E20) {
          distance = kInfLength;
          return;
        }
      }
      if (!inright) {
        boundary2 = unplaced.fRightVolume->DistanceToIn(hitpoint, dir);
        // Keep the same small progress convention for either constituent.
        d2 = Max(boundary2, Real_v(kHalfTolerance));
        if (d2 > 1E20) {
          distance = kInfLength;
          return;
        }
      }

      if (d1 > d2) {
        // propagate to left shape
        const Real_v boundaryDistance = snext + Max(boundary1, Real_v(0.));
        snext += d1;
        inleft = true;
        hitpoint += d1 * dir;

        // check if propagated point is inside right shape
        // check is done with a little push
        inright = unplaced.fRightVolume->Contains(hitpoint + kTolerance * dir);
        if (inright) {
          distance = boundaryDistance;
          return;
        }
        // here inleft=true, inright=false
      } else {
        // propagate to right shape
        const Real_v boundaryDistance = snext + Max(boundary2, Real_v(0.));
        snext += d2;
        inright = true;
        hitpoint += d2 * dir;

        // check if propagated point is inside left shape
        inleft = unplaced.fLeftVolume->Contains(hitpoint + kTolerance * dir);
        if (inleft) {
          distance = boundaryDistance;
          return;
        }
      }
      // here inleft=false, inright=true
    } // end while loop
    distance = snext;
    return;
  }

  // Collect both operands during their distance calls. Normal storage remains
  // local until the winner is known, so a losing operand cannot overwrite it.
  // Ordinary inlining breaks the forced-inline chain for nested Boolean hits.
  template <typename Real_v>
  VECCORE_ATT_HOST_DEVICE static void DistanceToOutWithHit(BooleanStruct const &unplaced, Vector3D<Real_v> const &point,
                                                           Vector3D<Real_v> const &direction, Real_v const &stepMax,
                                                           Real_v &distance, SurfaceHitView<Real_v> *hitInfo)
  {
    hitInfo->Clear();
    Vector3D<Real_v> leftNormal, rightNormal;
    SurfaceHitView<Real_v> leftHit{hitInfo->WantsNormal() ? &leftNormal : nullptr};
    SurfaceHitView<Real_v> rightHit{hitInfo->WantsNormal() ? &rightNormal : nullptr};
    const Real_v leftDistance  = BooleanHelper::DistanceToOutWithHit(unplaced.fLeftVolume, point, direction, stepMax,
                                                                     unplaced.fIdentityMask & 1, &leftHit);
    const Real_v rightDistance = BooleanHelper::DistanceToOutWithHit(unplaced.fRightVolume, point, direction, stepMax,
                                                                     unplaced.fIdentityMask & 2, &rightHit);
    distance                   = Min(leftDistance, rightDistance);

    // An intersection is a subset of either operand. At a tied exit, either
    // supporting incident plane suffices. Do not promote a negative invalid
    // operand answer (which may be a wrong-side sentinel) through the tie rule.
    const bool tied     = distance >= Real_v(0.) && Abs(leftDistance - rightDistance) <= Real_v(kTolerance);
    const bool useLeft  = leftHit.fLeftBehind && (leftDistance <= rightDistance || tied);
    const bool useRight = rightHit.fLeftBehind && (rightDistance <= leftDistance || tied);
    if (useLeft || useRight) {
      hitInfo->SetLeftBehind(true);
      if (hitInfo->WantsNormal()) hitInfo->SetNormal(useLeft ? leftNormal : rightNormal);
    }
  }

  /**
   * @brief Compute the distance from inside the intersection to leave it.
   * @param unplaced Boolean operands.
   * @param point Point in the Boolean frame.
   * @param direction Unit direction in the Boolean frame.
   * @param stepMax Distance limit forwarded to the operands.
   * @param[out] distance Minimum of the two operand exit distances.
   * @param[out] hitInfo Optional supporting-plane flag and Boolean-frame unit normal.
   * @details Leaving either constituent leaves the intersection. A winning
   * supporting operand normal also supports the intersection. Invalid results
   * leave normal storage untouched.
   */
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void DistanceToOut(BooleanStruct const &unplaced,
                                                                         Vector3D<Real_v> const &point,
                                                                         Vector3D<Real_v> const &direction,
                                                                         Real_v const &stepMax, Real_v &distance,
                                                                         SurfaceHitView<Real_v> *hitInfo = nullptr)
  {
    if (hitInfo) {
      DistanceToOutWithHit(unplaced, point, direction, stepMax, distance, hitInfo);
      return;
    }
    distance = Min(unplaced.fLeftVolume->PlacedDistanceToOut(point, direction, stepMax),
                   unplaced.fRightVolume->PlacedDistanceToOut(point, direction, stepMax));
  }

  /**
   * @brief Compute the safety from outside to enter the intersection.
   * @details If the point is already inside one constituent, only the other
   * constituent can block entry. Otherwise the closer constituent entry safety
   * is used, following the Geant4-style approximation kept by the legacy code.
   */
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SafetyToIn(BooleanStruct const &unplaced,
                                                                      Vector3D<Real_v> const &point, Real_v &safety)
  {
    // This is the Geant4 algorithm
    // TODO: ROOT seems to produce better safeties
    const auto insideA = unplaced.fLeftVolume->Contains(point);
    const auto insideB = unplaced.fRightVolume->Contains(point);

    if (!insideA && insideB) {
      safety = unplaced.fLeftVolume->SafetyToIn(point);
    } else {
      if (!insideB && insideA) {
        safety = unplaced.fRightVolume->SafetyToIn(point);
      } else {
        safety = Min(unplaced.fLeftVolume->SafetyToIn(point), unplaced.fRightVolume->SafetyToIn(point));
      }
    }
    return;
  }

  /**
   * @brief Compute the safety from inside the intersection to leave it.
   * @details Leaving either constituent leaves the intersection, so the smaller
   * constituent safety bounds the exit. A wrong-side constituent already returns
   * a non-positive safety; no separate classification is needed. Negative values
   * are clamped to zero, including on the wrong side. Builds that disable signed
   * constituent safeties retain explicit wrong-side classification.
   */
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SafetyToOut(BooleanStruct const &unplaced,
                                                                       Vector3D<Real_v> const &point, Real_v &safety)
  {
#ifdef VECGEOM_NO_WRONGSIDE_CONV
    // Some constituent kernels return unsigned distances in this configuration.
    if (unplaced.fLeftVolume->Inside(point) == kOutside || unplaced.fRightVolume->Inside(point) == kOutside) {
      safety = Real_v(0.);
      return;
    }
#endif
    // Placed-volume SafetyToOut expects constituent-local coordinates for both
    // operands; do not assume the left operand has identity placement.
    safety = unplaced.fLeftVolume->SafetyToOut(unplaced.fLeftVolume->GetTransformation()->Transform(point));
    // A non-positive left safety already fixes the clamped minimum at zero.
    if (safety > Real_v(0.))
      safety =
          Min(safety, unplaced.fRightVolume->SafetyToOut(unplaced.fRightVolume->GetTransformation()->Transform(point)));
    if (safety < Real_v(0.)) safety = Real_v(0.);
  }

  /**
   * @brief Compute an outward normal for the closest intersection surface.
   * @details The selected constituent is the one with the smaller safety to the
   * relevant boundary, with ties selecting the right constituent. Reuse the
   * local points returned by containment for the safety and normal queries,
   * then rotate the selected normal back to the intersection frame.
   */
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void NormalKernel(BooleanStruct const &unplaced,
                                                                        Vector3D<Real_v> const &point,
                                                                        Vector3D<Real_v> &normal, bool &valid)
  {
    valid = false;

    VPlacedVolume const *const fPtrSolidA = unplaced.fLeftVolume;
    VPlacedVolume const *const fPtrSolidB = unplaced.fRightVolume;
    Vector3D<Real_v> pointA, pointB;
    Real_v safetyA, safetyB;

    if (fPtrSolidA->Contains(point, pointA)) {
      // Placed-volume SafetyToOut already expects constituent-local coordinates.
      safetyA = fPtrSolidA->SafetyToOut(pointA);
    } else {
      safetyA = fPtrSolidA->GetUnplacedVolume()->SafetyToIn(pointA);
    }

    if (fPtrSolidB->Contains(point, pointB)) {
      safetyB = fPtrSolidB->SafetyToOut(pointB);
    } else {
      safetyB = fPtrSolidB->GetUnplacedVolume()->SafetyToIn(pointB);
    }
    const bool onA          = safetyA < safetyB;
    auto const *const owner = onA ? fPtrSolidA : fPtrSolidB;
    Vector3D<Real_v> localNormal;
    valid = owner->GetUnplacedVolume()->Normal(onA ? pointA : pointB, localNormal);
    owner->GetTransformation()->InverseTransformDirection(localNormal, normal);
  }
}; // End struct BooleanImplementation

} // namespace VECGEOM_IMPL_NAMESPACE

} // namespace vecgeom

#endif /* BooleanImplementation_H_ */
