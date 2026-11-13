/*
 * PolyconeImplementation.h
 *
 *  Created on: Dec 8, 2014
 *      Author: swenzel
 */

/// @file PolyconeImplementation.h
/// @brief Geometry kernels for the classic polycone solid.

/// History notes:
/// Jan-March 2017: revision + moving to use new Cone Kernels (Raman Sehgal)
/// May-June 2017: revision + moving to new Structure (Raman Sehgal)

#ifndef VECGEOM_VOLUMES_KERNEL_POLYCONEIMPLEMENTATION_H_
#define VECGEOM_VOLUMES_KERNEL_POLYCONEIMPLEMENTATION_H_

#include "VecGeom/base/Vector3D.h"
#include "VecGeom/volumes/PolyconeStruct.h"
#include "VecGeom/volumes/SurfaceHitView.h"
#include "VecGeom/volumes/kernel/GenericKernels.h"
#include <VecCore/VecCore>
#include "VecGeom/volumes/kernel/ConeImplementation.h"
#include "VecGeom/volumes/kernel/shapetypes/ConeTypes.h"

#include <cstdio>

namespace vecgeom {

VECGEOM_DEVICE_DECLARE_CONV_TEMPLATE(struct, PolyconeImplementation, typename);

inline namespace VECGEOM_IMPL_NAMESPACE {

template <typename T>
class SPlacedPolycone;
template <typename T>
class SUnplacedPolycone;

template <typename polyconeTypeT>
/**
 * @brief Kernel entry points for the classic polycone.
 *
 * The classic polycone is implemented as a union of shifted cone sections plus
 * an optional phi wedge. Most point and distance queries are delegated to the
 * section cones, while shared-plane hand-off and contour-based safety are
 * resolved at the polycone layer.
 */
struct PolyconeImplementation {

  using UnplacedStruct_t = PolyconeStruct<Precision>;
  using UnplacedVolume_t = SUnplacedPolycone<polyconeTypeT>;
  using PlacedShape_t    = SPlacedPolycone<UnplacedVolume_t>;
  using SectionStruct_t  = RadialConeStruct<Precision>;
  using SectionCone_t    = ConeTypes::UniversalConeWithoutPhi;
  using SectionImpl_t    = ConeImplementation<SectionCone_t, SectionStruct_t>;

  /** Record a phi-plane exit, which supports the solid only for an opening no larger than pi. */
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SetPhiExitHit(UnplacedStruct_t const &polycone,
                                                                         bool startSurface,
                                                                         SurfaceHitView<Real_v> *hitInfo)
  {
    const bool leftBehind = polycone.fDeltaPhi <= kPi;
    hitInfo->SetLeftBehind(leftBehind);
    if (leftBehind && hitInfo->WantsNormal()) {
      const Vector3D<Precision> normal =
          startSurface ? polycone.fPhiWedge.GetNormal1() : polycone.fPhiWedge.GetNormal2();
      hitInfo->SetNormal(Vector3D<Real_v>(-normal.x(), -normal.y(), Real_v(0.)));
    }
  }

  /** Record a selected outer-section exit using its cached whole-contour support classification. */
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SetOuterExitHit(PolyconeSection const &section,
                                                                           Vector3D<Real_v> const &hitPoint,
                                                                           SurfaceHitView<Real_v> *hitInfo)
  {
    hitInfo->SetLeftBehind(section.fConvex);
    if (section.fConvex && hitInfo->WantsNormal()) {
      const Vector3D<Real_v> localPoint = hitPoint - Vector3D<Real_v>(Real_v(0.), Real_v(0.), section.fShift);
      // Reuse the section's cached unit-normal factors; only the radial direction needs normalization.
      const Real_v invRhoSec = section.fSolid.fInvSecRMax / NonZero(localPoint.Perp());
      hitInfo->SetNormal(
          Vector3D<Real_v>(localPoint.x() * invRhoSec, localPoint.y() * invRhoSec, section.fSolid.fZNormOuter));
    }
  }

  /** Return whether the final hit point is incident on the selected section's outer conical face. */
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static bool IsIncidentOnOuterSurface(PolyconeSection const &section,
                                                                                    Vector3D<Real_v> const &hitPoint)
  {
    const Vector3D<Real_v> localPoint = hitPoint - Vector3D<Real_v>(Real_v(0.), Real_v(0.), section.fShift);
    const Real_v outerRadius = ConeUtilities::GetRadiusOfConeAtPoint<Real_v, false>(section.fSolid, localPoint.z());
    return localPoint.Perp2() >= MakeMinusTolerantSquare<true>(outerRadius, section.fSolid.fOuterTolerance);
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void AddPhiNormal(
      UnplacedStruct_t const &polycone, PolyconeSection const &section, Vector3D<Real_v> const &point,
      bool const onStart, bool const onEnd, Vector3D<Real_v> &normal, int &surfaceCount, bool &valid)
  {
    if (!(onStart || onEnd)) return;

    bool radialContains = false;
    SectionImpl_t::Contains(section.fSolid, point - Vector3D<Real_v>(Real_v(0.), Real_v(0.), Real_v(section.fShift)),
                            radialContains);
    if (!radialContains) return;

    if (surfaceCount == 0) normal.Set(Real_v(0.), Real_v(0.), Real_v(0.));
    if (onStart) {
      normal -= Vector3D<Real_v>(polycone.fPhiWedge.GetNormal1());
      ++surfaceCount;
    }
    if (onEnd) {
      normal -= Vector3D<Real_v>(polycone.fPhiWedge.GetNormal2());
      ++surfaceCount;
    }
    valid = surfaceCount <= 2;
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static bool SurfaceNormal(UnplacedStruct_t const &polycone,
                                                                         Vector3D<Real_v> const &point,
                                                                         Vector3D<Real_v> &normal)
  {
    if (polycone.fDeltaPhi < kTwoPi && !polycone.fPhiWedge.ContainsWithBoundary(point)) return false;

    int index = polycone.GetSectionIndex(point.z() - Real_v(kTolerance));
    if (index < 0) {
      if (index == -1) normal.Set(Real_v(0.), Real_v(0.), Real_v(-1.));
      if (index == -2) normal.Set(Real_v(0.), Real_v(0.), Real_v(1.));
      return true;
    }

    const PolyconeSection &section = polycone.GetSection(index);
    Vector3D<Precision> radialNormal;
    int surfaceCount = 0;
    bool valid       = section.fSolid.NormalWithCount(
        Vector3D<Precision>(point) - Vector3D<Precision>(0., 0., section.fShift), radialNormal, surfaceCount);
    normal = Vector3D<Real_v>(radialNormal);
    const bool onStart =
        polycone.fDeltaPhi < kTwoPi && polycone.fPhiWedge.template IsOnSurfaceGeneric<Real_v, true>(point);
    const bool onEnd =
        polycone.fDeltaPhi < kTwoPi && polycone.fPhiWedge.template IsOnSurfaceGeneric<Real_v, false>(point);
    AddPhiNormal(polycone, section, point, onStart, onEnd, normal, surfaceCount, valid);
    if (surfaceCount > 1) normal = normal.Unit();

    if (index + 1 < polycone.GetNSections() &&
        Abs(point.z() - Real_v(polycone.GetZAtPlane(index + 1))) < Real_v(kTolerance)) {
      const PolyconeSection &nextSection = polycone.GetSection(index + 1);
      Vector3D<Precision> nextNormal;
      int nextSurfaceCount = 0;
      const bool nextValid = nextSection.fSolid.NormalWithCount(
          Vector3D<Precision>(point) - Vector3D<Precision>(0., 0., nextSection.fShift), nextNormal, nextSurfaceCount);
      Vector3D<Real_v> nextSectionNormal(nextNormal);
      bool nextSectionValid = nextValid;
      AddPhiNormal(polycone, nextSection, point, onStart, onEnd, nextSectionNormal, nextSurfaceCount, nextSectionValid);
      if (nextSurfaceCount > 1) nextSectionNormal = nextSectionNormal.Unit();
      const bool hasNormal     = surfaceCount > 0;
      const bool nextHasNormal = nextSurfaceCount > 0;
      if (!hasNormal && nextHasNormal) {
        normal = nextSectionNormal;
        valid  = nextSectionValid;
      } else if (hasNormal && nextHasNormal) {
        bool nextContains = false;
        SectionImpl_t::Contains(nextSection.fSolid,
                                point + Real_v(10. * kTolerance) * Vector3D<Real_v>(nextNormal) -
                                    Vector3D<Real_v>(0., 0., Real_v(nextSection.fShift)),
                                nextContains);
        if (nextContains) {
          normal = nextSectionNormal;
          valid  = nextSectionValid;
        } else {
          normal += nextSectionNormal;
          valid = valid && nextSectionValid;
          if (normal.Mag2() < Real_v(kTolerance)) normal.Set(point.x(), point.y(), Real_v(0.));
        }
      }
    }

    if (normal.Mag2() > Real_v(0.)) normal = normal.Unit();
    return valid;
  }

  /**
   * @brief Return the squared minimum `z` separation between a query `z` and an
   * `(r,z)` segment.
   *
   * This is the cheap directional pruning bound used by the contour walk in
   * `SafetyToIn`.
   */
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static Real_v SegmentMinDeltaZSquared(Real_v const &z,
                                                                                     Vector3D<Precision> const &a,
                                                                                     Vector3D<Precision> const &b)
  {
    const Real_v z0 = Min(Real_v(a.y()), Real_v(b.y()));
    const Real_v z1 = Max(Real_v(a.y()), Real_v(b.y()));
    if (z < z0) {
      const Real_v dz = z0 - z;
      return dz * dz;
    }
    if (z > z1) {
      const Real_v dz = z - z1;
      return dz * dz;
    }
    return Real_v(0.);
  }

  /**
   * @brief Return the squared distance from `(rho,z)` to one `(r,z)` segment.
   *
   * The projection is clamped to the segment end points. The result stays in
   * squared form so the contour search can avoid intermediate square roots.
   */
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static Real_v DistanceSquaredToRZSegment(Real_v const &rho,
                                                                                        Real_v const &z,
                                                                                        Vector3D<Precision> const &a,
                                                                                        Vector3D<Precision> const &b)
  {
    const Real_v ax = Real_v(a.x());
    const Real_v az = Real_v(a.y());
    const Real_v bx = Real_v(b.x());
    const Real_v bz = Real_v(b.y());
    const Real_v dx = bx - ax;
    const Real_v dz = bz - az;
    const Real_v l2 = dx * dx + dz * dz;

    if (l2 <= Real_v(0.)) {
      const Real_v qx = rho - ax;
      const Real_v qz = z - az;
      return qx * qx + qz * qz;
    }

    Real_v t = ((rho - ax) * dx + (z - az) * dz) / l2;
    t        = Max(Real_v(0.), Min(Real_v(1.), t));

    const Real_v rx = rho - (ax + t * dx);
    const Real_v rz = z - (az + t * dz);
    return rx * rx + rz * rz;
  }

  /**
   * @brief Evaluate one cone section for containment or inside classification.
   *
   * @param unplaced Polycone data model.
   * @param isect Section index. Negative indices are treated as fully outside.
   * @param polyconePoint Query point in polycone coordinates.
   * @param secFullyInside Output flag for section containment.
   * @param secFullyOutside Output flag for section exclusion.
   */
  template <typename Real_v, bool ForInside>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void GenericKernelForASection(
      UnplacedStruct_t const &unplaced, int isect, Vector3D<Real_v> const &polyconePoint, bool &secFullyInside,
      bool &secFullyOutside)
  {

    // using namespace PolyconeTypes;
    using namespace ConeTypes;

    if (isect < 0) {
      secFullyInside  = false;
      secFullyOutside = true;
      return;
    }

    PolyconeSection const &sec    = unplaced.GetSection(isect);
    Vector3D<Precision> secLocalp = polyconePoint - Vector3D<Precision>(0, 0, sec.fShift);

    ConeHelpers<Real_v, SectionCone_t, SectionStruct_t>::template GenericKernelForContainsAndInside<ForInside>(
        sec.fSolid, secLocalp, secFullyInside, secFullyOutside);
  }

  template <typename Real_v, bool ForInside>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static bool ClassifyOrdinarySection(UnplacedStruct_t const &polycone,
                                                                                   Vector3D<Real_v> const &point,
                                                                                   bool &completelyInside,
                                                                                   bool &completelyOutside)
  {
    const int section = polycone.GetSectionIndex(point.z());
    if (section < 0 || Abs(point.z() - Real_v(polycone.GetZAtPlane(section))) <= Real_v(kTolerance) ||
        Abs(point.z() - Real_v(polycone.GetZAtPlane(section + 1))) <= Real_v(kTolerance))
      return false;
    GenericKernelForASection<Real_v, ForInside>(polycone, section, point, completelyInside, completelyOutside);
    return true;
  }

  /**
   * @brief Check whether a point belongs to the polycone.
   *
   * This is the boundary-inclusive predicate used by `Contains`.
   */
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void Contains(UnplacedStruct_t const &polycone,
                                                                    Vector3D<Real_v> const &point, bool &inside)
  {
    const int sectionCount = polycone.GetNSections();
    if (sectionCount == 1) {
      bool unused = false, outside = false;
      GenericKernelForASection<Real_v, false>(polycone, 0, point, unused, outside);
      if (!outside && polycone.fDeltaPhi < kTwoPi) {
        bool wedgeInside = false, wedgeOutside = false;
        polycone.fPhiWedge.template GenericKernelForContainsAndInside<false>(point, wedgeInside, wedgeOutside);
        outside |= wedgeOutside;
      }
      inside = !outside;
      return;
    }
    bool unused  = false;
    bool outside = false;
    if (ClassifyOrdinarySection<Real_v, false>(polycone, point, unused, outside)) {
      if (!outside && polycone.fDeltaPhi < kTwoPi) {
        bool wedgeInside = false, wedgeOutside = false;
        polycone.fPhiWedge.template GenericKernelForContainsAndInside<false>(point, wedgeInside, wedgeOutside);
        outside |= wedgeOutside;
      }
      inside = !outside;
      return;
    }
    GenericKernelForContainsAndInside<Real_v, false>(polycone, point, unused, outside);
    if (!outside && polycone.fDeltaPhi < kTwoPi) {
      bool wedgeInside = false, wedgeOutside = false;
      polycone.fPhiWedge.template GenericKernelForContainsAndInside<false>(point, wedgeInside, wedgeOutside);
      outside |= wedgeOutside;
    }
    inside = !outside;
  }

  /**
   * @brief Classify a point as `kInside`, `kSurface`, or `kOutside`.
   *
   * Shared section planes are resolved at the polycone level before mapping the
   * result to `EInside`.
   */
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void Inside(UnplacedStruct_t const &polycone,
                                                                  Vector3D<Real_v> const &point, Inside_t &inside)
  {
    // A very narrow wedge rejects enough points to offset its more expensive
    // surface classification. Wider sectors keep the cheaper radial rejection
    // first; the shape-level branch is uniform for every query.
    constexpr Precision kPhiFirstThreshold = 0.4;
    bool wedgeInside = true, wedgeOutside = false;
    const bool phiFirst = polycone.fDeltaPhi < kPhiFirstThreshold;
    if (phiFirst) {
      polycone.fPhiWedge.template GenericKernelForContainsAndInside<true>(point, wedgeInside, wedgeOutside);
      if (wedgeOutside) {
        inside = EInside::kOutside;
        return;
      }
    }

    const int sectionCount = polycone.GetNSections();
    if (sectionCount == 1) {
      bool completelyinside = false, completelyoutside = false;
      GenericKernelForASection<Real_v, true>(polycone, 0, point, completelyinside, completelyoutside);
      if (!completelyoutside && polycone.fDeltaPhi < kTwoPi) {
        if (!phiFirst)
          polycone.fPhiWedge.template GenericKernelForContainsAndInside<true>(point, wedgeInside, wedgeOutside);
        completelyoutside |= wedgeOutside;
        completelyinside &= wedgeInside;
      }
      inside = completelyoutside ? EInside::kOutside : (completelyinside ? EInside::kInside : EInside::kSurface);
      return;
    }

    bool completelyinside  = false;
    bool completelyoutside = false;
    if (ClassifyOrdinarySection<Real_v, true>(polycone, point, completelyinside, completelyoutside)) {
      if (!completelyoutside && polycone.fDeltaPhi < kTwoPi) {
        if (!phiFirst)
          polycone.fPhiWedge.template GenericKernelForContainsAndInside<true>(point, wedgeInside, wedgeOutside);
        completelyoutside |= wedgeOutside;
        completelyinside &= wedgeInside;
      }
      inside = completelyoutside ? EInside::kOutside : (completelyinside ? EInside::kInside : EInside::kSurface);
      return;
    }
    GenericKernelForContainsAndInside<Real_v, true>(polycone, point, completelyinside, completelyoutside);
    if (!completelyoutside && polycone.fDeltaPhi < kTwoPi) {
      if (!phiFirst)
        polycone.fPhiWedge.template GenericKernelForContainsAndInside<true>(point, wedgeInside, wedgeOutside);
      completelyoutside |= wedgeOutside;
      completelyinside &= wedgeInside;
    }
    inside = EInside::kSurface;
    if (completelyoutside) inside = EInside::kOutside;
    if (completelyinside) inside = EInside::kInside;
  }

  /**
   * @brief Shared kernel for `Contains` and `Inside`.
   *
   * The point is matched against the section owning `z - kTolerance` and
   * `z + kTolerance`. This lets the polycone distinguish ordinary section
   * ownership from shared-plane cases near repeated-`z` boundaries.
   */
  template <typename Real_v, bool ForInside>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void GenericKernelForContainsAndInside(
      UnplacedStruct_t const &unplaced, Vector3D<Real_v> const &localPoint, bool &completelyInside,
      bool &completelyOutside)
  {
    int indexLow  = unplaced.GetSectionIndex(localPoint.z() - kTolerance);
    int indexHigh = unplaced.GetSectionIndex(localPoint.z() + kTolerance);
    if (indexLow < 0 && indexHigh < 0) {
      completelyOutside = true;
      return;
    }
    if (indexLow < 0 && indexHigh == 0) {
      // Check location in section 0 and return
      GenericKernelForASection<Real_v, ForInside>(unplaced, 0, localPoint, completelyInside, completelyOutside);
      return;
    }
    if (indexHigh < 0 && indexLow == (unplaced.GetNSections() - 1)) {
      // Check location in section N-1 and return
      GenericKernelForASection<Real_v, ForInside>(unplaced, (unplaced.GetNSections() - 1), localPoint, completelyInside,
                                                  completelyOutside);

      return;
    }
    if (indexLow >= 0 && indexHigh >= 0) {
      if (indexLow == indexHigh) {
        // Check location in section indexLow and return
        GenericKernelForASection<Real_v, ForInside>(unplaced, indexLow, localPoint, completelyInside,
                                                    completelyOutside);

        return;
      } else {

        bool secInLow = false, secOutLow = false;
        bool secInHigh = false, secOutHigh = false;
        GenericKernelForASection<Real_v, ForInside>(unplaced, indexLow, localPoint, secInLow, secOutLow);
        GenericKernelForASection<Real_v, ForInside>(unplaced, indexHigh, localPoint, secInHigh, secOutHigh);
        bool surfLow  = !secInLow && !secOutLow;
        bool surfHigh = !secInHigh && !secOutHigh;

        if (surfLow && surfHigh) {
          if constexpr (ForInside) {
            const bool radialLow  = SectionPointOnConicalSurface(unplaced.GetSection(indexLow), localPoint);
            const bool radialHigh = SectionPointOnConicalSurface(unplaced.GetSection(indexHigh), localPoint);
            if (radialLow || radialHigh) {
              // Adjacent sections can both see a point as surface near a
              // shared z plane. Keep exposed conical edge/side hits as
              // kSurface; only pure internal cap overlap collapses to inside.
              return;
            }
          }
          // A point can be surface for both adjacent sections on a shared
          // plane. Containment keeps the cheap "inside the union" answer here;
          // entry/exit hand-off is refined later by the distance kernels using
          // the two section end annuli.
          completelyInside = true;
          return;
        } else {
          // else if point is on surface of only one of the two sections then point is actually on surface , the default
          // case,
          // so no need to check

          // What needs to check is if it is outside both ie. Outside indexLow section and Outside indexHigh section
          // then it is certainly outside
          if (secOutLow && secOutHigh) {
            completelyOutside = true;
            return;
          }
        }
      }
    }
  }

  /**
   * @brief Return whether a point lies on one section's conical radial surface.
   * @details This is used only in rare shared-plane classification, where two
   * adjacent sections both report surface and the polycone must distinguish an
   * internal cap overlap from an exposed conical edge/side hit.
   */
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static bool SectionPointOnConicalSurface(PolyconeSection const &section,
                                                                                        Vector3D<Real_v> const &point)
  {
    const Vector3D<Real_v> localPoint = point - Vector3D<Precision>(0, 0, section.fShift);
    if (Abs(SectionImpl_t::template SafeDistanceToConicalSurface<Real_v, false>(section.fSolid, localPoint)) <
        Real_v(section.fSolid.fOuterTolerance)) {
      return true;
    }
    if (ConeTypes::checkRminTreatment<polyconeTypeT>(section.fSolid)) {
      return Abs(SectionImpl_t::template SafeDistanceToConicalSurface<Real_v, true>(section.fSolid, localPoint)) <
             Real_v(section.fSolid.fInnerTolerance);
    }
    return false;
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SectionDistanceToIn(PolyconeSection const &section,
                                                                               evolution::Wedge const *wedge,
                                                                               Vector3D<Real_v> const &point,
                                                                               Vector3D<Real_v> const &direction,
                                                                               Real_v const &stepMax, Real_v &distance)
  {
    SectionImpl_t::template DistanceToIn<Real_v>(section.fSolid, point, direction, stepMax, distance, wedge);
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void DistanceToIn(UnplacedStruct_t const &polycone,
                                                                        Vector3D<Real_v> const &point,
                                                                        Vector3D<Real_v> const &direction,
                                                                        Real_v const &stepMax, Real_v &distance)
  {
    // Beyond the outer radial envelope, an outward ray cannot reach any section.
    if (polycone.GetNSections() > 1 && point.Perp2() > Real_v(polycone.fOuterBoundSquared) &&
        point.x() * direction.x() + point.y() * direction.y() >= Real_v(0.)) {
      distance = kInfLength;
      return;
    }
    if (polycone.fDeltaPhi >= kTwoPi) {
      RadialDistanceToIn(polycone, point, direction, stepMax, distance);
      return;
    }
    RadialDistanceToIn(polycone, point, direction, stepMax, distance, &polycone.fPhiWedge);
  }

  template <typename Real_v>
#ifdef VECCORE_CUDA_DEVICE_COMPILATION
  VECCORE_FORCE_NOINLINE
#else
  VECGEOM_FORCE_INLINE
#endif
      VECCORE_ATT_HOST_DEVICE static void RadialDistanceToIn(UnplacedStruct_t const &polycone,
                                                             Vector3D<Real_v> const &point,
                                                             Vector3D<Real_v> const &direction, Real_v const &stepMax,
                                                             Real_v &distance, evolution::Wedge const *wedge = nullptr)
  {
    // using namespace PolyconeTypes;
    Vector3D<Real_v> p = point;
    Vector3D<Real_v> v = direction;

    // TODO: add bounding box check maybe??

    distance      = kInfLength;
    int increment = (v.z() > 0) ? 1 : -1;
    if (Abs(v.z()) < Real_v(kTolerance)) increment = 0;
    int index = polycone.GetSectionIndex(p.z());
    if (index == -1) index = 0;
    if (index == -2) index = polycone.GetNSections() - 1;

    if (increment) {
      // Keep the rare z-grazing hand-off path out of the common section walk.
      bool rejectedHandoff = false;
      do {
        const bool rejectCurrentEntryCap = rejectedHandoff;
        rejectedHandoff                  = false;

        // now we have to find a section
        PolyconeSection const &sec = polycone.GetSection(index);

        SectionDistanceToIn(sec, wedge, p - Vector3D<Precision>(0, 0, sec.fShift), v, stepMax, distance);

        if (distance < kInfLength) {
          const int nextIndex     = index + increment;
          const Real_v localHitZ  = p.z() + distance * v.z() - Real_v(sec.fShift);
          const Real_v touchedEnd = increment > 0 ? Real_v(sec.fSolid.fDz) : Real_v(-sec.fSolid.fDz);
          if (rejectCurrentEntryCap) {
            // The previous section touched this section's entry cap without a
            // material interval ahead. Do not let the same cap hit be accepted
            // again by the next section cone.
            const Real_v entryEnd = -touchedEnd;
            if (Abs(localHitZ - entryEnd) < Real_v(kTolerance)) {
              const Vector3D<Real_v> hit = p + distance * v;
              const Real_v rho2          = hit.x() * hit.x() + hit.y() * hit.y();
              if (!SectionEndHasForwardMaterial(sec, hit, rho2, increment < 0, v)) distance = kInfLength;
            }
          }
          if (distance < kInfLength && Abs(localHitZ - touchedEnd) < Real_v(kTolerance)) {
            if (nextIndex >= 0 && nextIndex < polycone.GetNSections()) {
              const Vector3D<Real_v> hit = p + distance * v;
              const Real_v rho2          = hit.x() * hit.x() + hit.y() * hit.y();

              const PolyconeSection &nextSec = polycone.GetSection(nextIndex);
              const Precision nextRMin       = increment > 0 ? nextSec.fSolid.fRmin1 : nextSec.fSolid.fRmin2;
              const Precision nextRMax       = increment > 0 ? nextSec.fSolid.fRmax1 : nextSec.fSolid.fRmax2;

              const Real_v nextMin   = Max(Real_v(0.), Real_v(nextRMin) - Real_v(kTolerance));
              const Real_v nextMax   = Real_v(nextRMax) + Real_v(kTolerance);
              const Real_v nextMinSq = nextMin * nextMin;
              const Real_v nextMaxSq = nextMax * nextMax;

              const bool onNextEnd = (rho2 > nextMinSq) && (rho2 < nextMaxSq);
              const bool ownsForwardInterval =
                  onNextEnd && SectionEndHasForwardMaterial(nextSec, hit, rho2, increment < 0, v);

              if (!ownsForwardInterval) {
                rejectedHandoff = onNextEnd;
                distance        = kInfLength;
              } else {
                // A section-end hit from DistanceToIn is an edge touch of the
                // current cone. Accept it only when the next section owns a
                // non-zero material interval along the ray; otherwise ignore
                // this grazing edge and continue walking sections.
                return;
              }
            } else {
              distance = kInfLength;
            }
          }
        }

        if (distance < kInfLength) break;
        index += increment;
      } while (index >= 0 && index < polycone.GetNSections());
      return;
    }

    // No z-section increment means this is a z-grazing ray: it is parallel
    // or nearly parallel to section end planes.
    const bool onBottomCap = Abs(p.z() - Real_v(polycone.GetZAtPlane(0))) < Real_v(kTolerance);
    const bool onTopCap    = Abs(p.z() - Real_v(polycone.GetZAtPlane(polycone.GetNSections()))) < Real_v(kTolerance);
    if ((onBottomCap && v.z() <= Real_v(0.)) || (onTopCap && v.z() >= Real_v(0.))) {
      // Near-horizontal rays moving out of a real terminal cap, including
      // exact grazes, do not own a forward material interval for
      // DistanceToIn. Inward near-grazing cap rays are still handled by the
      // section cone and must be paired with a non-zero DistanceToOut
      // continuation.
      return;
    }

    PolyconeSection const &sec = polycone.GetSection(index);
    SectionDistanceToIn(sec, wedge, p - Vector3D<Precision>(0, 0, sec.fShift), v, stepMax, distance);

    const bool zPositive = v.z() > Real_v(0.);
    const bool zNegative = v.z() < Real_v(0.);
    if (index + 1 < polycone.GetNSections() &&
        Abs(p.z() - Real_v(polycone.GetZAtPlane(index + 1))) < Real_v(kTolerance)) {
      Real_v nextDistance            = kInfLength;
      const PolyconeSection &nextSec = polycone.GetSection(index + 1);
      // On a shared z plane, z-grazing rays cannot pick a reliable section
      // from the section index alone. Positive-z rays use the high-z
      // section; exactly horizontal rays use the farther finite circle
      // because the nearer one can be only a grazing boundary of the
      // exposed transition annulus.
      SectionDistanceToIn(nextSec, wedge, p - Vector3D<Precision>(0, 0, nextSec.fShift), v, stepMax, nextDistance);
      if (zPositive) {
        distance = nextDistance;
      } else if (!zNegative && nextDistance < kInfLength) {
        distance = (distance < kInfLength) ? Max(distance, nextDistance) : nextDistance;
      }
    }
    if (index > 0 && Abs(p.z() - Real_v(polycone.GetZAtPlane(index))) < Real_v(kTolerance)) {
      Real_v previousDistance            = kInfLength;
      const PolyconeSection &previousSec = polycone.GetSection(index - 1);
      // Symmetric fallback for a point selected on the high-z side of a
      // shared plane: negative-z rays use the low-z section, while exactly
      // horizontal rays keep the material-after, farther finite entry.
      SectionDistanceToIn(previousSec, wedge, p - Vector3D<Precision>(0, 0, previousSec.fShift), v, stepMax,
                          previousDistance);
      if (zNegative) {
        distance = previousDistance;
      } else if (!zPositive && previousDistance < kInfLength) {
        distance = (distance < kInfLength) ? Max(distance, previousDistance) : previousDistance;
      }
    }
    return;
  }

  /**
   * @brief Return whether one section end annulus owns a point.
   * @details `rho2` is passed by the caller so shared-plane resolution can test
   * both neighboring sections without recomputing the same radial coordinate.
   * Ownership is a local endpoint predicate: the point must lie inside the
   * end-cap radial band and, for phi-limited sections, inside the boundary-
   * inclusive phi wedge.
   */
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static bool SectionEndOwnsPoint(PolyconeSection const &section,
                                                                               Vector3D<Real_v> const &point,
                                                                               Real_v const &rho2, bool upperEnd)
  {
    const Precision rmin = upperEnd ? section.fSolid.fRmin2 : section.fSolid.fRmin1;
    const Precision rmax = upperEnd ? section.fSolid.fRmax2 : section.fSolid.fRmax1;

    const Real_v minRadius = Max(Real_v(0.), Real_v(rmin) - Real_v(kTolerance));
    const Real_v maxRadius = Real_v(rmax) + Real_v(kTolerance);
    const bool insideInner = (rmin <= kTolerance) || (rho2 > minRadius * minRadius);
    const bool insideOuter = rho2 < maxRadius * maxRadius;
    if (!(insideInner && insideOuter)) return false;

    return true;
  }

  /**
   * @brief Return whether an owned section end has material immediately ahead.
   * @details The caller must already have established endpoint ownership. The
   * derivative checks reject edge touches that remain outside the conical
   * radial interval after an infinitesimal forward step, without paying for a
   * full section `Inside` query.
   */
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static bool SectionEndHasForwardMaterial(PolyconeSection const &section,
                                                                                        Vector3D<Real_v> const &point,
                                                                                        Real_v const &rho2,
                                                                                        bool upperEnd,
                                                                                        Vector3D<Real_v> const &dir)
  {
    const bool entersFromLower = !upperEnd && dir.z() > Real_v(kTolerance);
    const bool entersFromUpper = upperEnd && dir.z() < -Real_v(kTolerance);
    if (!(entersFromLower || entersFromUpper)) return false;

    const Real_v rmin = upperEnd ? Real_v(section.fSolid.fRmin2) : Real_v(section.fSolid.fRmin1);
    const Real_v rmax = upperEnd ? Real_v(section.fSolid.fRmax2) : Real_v(section.fSolid.fRmax1);
    bool nearInner    = false;

    if (ConeTypes::checkRminTreatment<polyconeTypeT>(section.fSolid)) {
      const Real_v innerLimit = rmin + Real_v(kTolerance);
      nearInner               = rho2 < innerLimit * innerLimit;
    }

    const Real_v outerLimit = Max(Real_v(0.), rmax - Real_v(kTolerance));
    const bool nearOuter    = rho2 > outerLimit * outerLimit;
    if (!(nearInner || nearOuter)) return true;

    const Real_v radialDerivative = Real_v(2.) * (point.x() * dir.x() + point.y() * dir.y());
    const Real_v zDerivative      = dir.z();

    if (nearInner) {
      const Real_v innerDerivative = Real_v(2.) * rmin * Real_v(section.fSolid.fInnerSlope) * zDerivative;
      if (radialDerivative <= innerDerivative + Real_v(kHalfTolerance)) return false;
    }

    if (nearOuter) {
      const Real_v outerDerivative = Real_v(2.) * rmax * Real_v(section.fSolid.fOuterSlope) * zDerivative;
      if (radialDerivative >= outerDerivative - Real_v(kHalfTolerance)) return false;
    }

    return true;
  }

  /**
   * @brief Return whether a section owns material immediately after an end plane.
   * @details Endpoint ownership alone is not enough when an inner radius opens
   * or an outer radius closes at the end plane: the endpoint may be on the
   * boundary, while the forward ray immediately leaves the section material.
   */
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static bool SectionEndOwnsForwardInterval(PolyconeSection const &section,
                                                                                         Vector3D<Real_v> const &point,
                                                                                         Real_v const &rho2,
                                                                                         bool upperEnd,
                                                                                         Vector3D<Real_v> const &dir)
  {
    if (!SectionEndOwnsPoint(section, point, rho2, upperEnd)) return false;
    return SectionEndHasForwardMaterial(section, point, rho2, upperEnd, dir);
  }

  /**
   * @brief Select the section owning a point on a shared z plane.
   * @details The low section is tested on its upper end and the high section on
   * its lower end. Unique ownership handles exposed transition annuli. If both
   * sections own the point, the z direction selects the immediate continuation;
   * exactly z-grazing rays keep the low section, matching the legacy section
   * ordering convention.
   */
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static bool ResolveSharedZSection(UnplacedStruct_t const &polycone,
                                                                                 Vector3D<Real_v> const &point,
                                                                                 Vector3D<Real_v> const &dir,
                                                                                 Real_v const &rho2, int lowIndex,
                                                                                 int highIndex, int &index)
  {
    const bool lowOwns  = SectionEndOwnsPoint(polycone.GetSection(lowIndex), point, rho2, true);
    const bool highOwns = SectionEndOwnsPoint(polycone.GetSection(highIndex), point, rho2, false);

    if (lowOwns && !highOwns) {
      index = lowIndex;
      return true;
    }
    if (!lowOwns && highOwns) {
      index = highIndex;
      return true;
    }
    if (!lowOwns && !highOwns) {
      const bool useLowConical  = SectionPointOnConicalSurface(polycone.GetSection(lowIndex), point);
      const bool useHighConical = SectionPointOnConicalSurface(polycone.GetSection(highIndex), point);
      if (useLowConical || useHighConical) {
        // Steep near-repeated sections can put a real conical surface hit just
        // outside the endpoint-radius band while still inside the z tolerance
        // of the shared plane. Let the conical side own such starts; the
        // section cone decides whether this is a zero exit or a finite
        // continuation.
        index = useLowConical ? lowIndex : highIndex;
        return true;
      }
      return false;
    }

    // When both sections own the shared plane, prefer the section that owns a
    // forward material interval. Otherwise keep the section being left so the
    // distance query returns the immediate boundary exit.
    if (dir.z() > Real_v(kTolerance)) {
      index =
          SectionEndHasForwardMaterial(polycone.GetSection(highIndex), point, rho2, false, dir) ? highIndex : lowIndex;
    } else if (dir.z() < -Real_v(kTolerance)) {
      index =
          SectionEndHasForwardMaterial(polycone.GetSection(lowIndex), point, rho2, true, dir) ? lowIndex : highIndex;
    } else {
      index = lowIndex;
    }
    return true;
  }

  /**
   * @brief Compute the distance from an inside point to the next polycone exit.
   *
   * The starting section is selected from the local point and direction. The
   * query is then delegated to section cones while the polycone layer manages
   * shared-plane transitions between neighboring sections. When requested,
   * hit information is classified only after the final section is known.
   */
  template <typename Real_v, bool RecordHitInfo>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void DistanceToOutImpl(UnplacedStruct_t const &polycone,
                                                                             Vector3D<Real_v> const &point,
                                                                             Vector3D<Real_v> const &dir,
                                                                             Real_v const &stepMax, Real_v &distance,
                                                                             SurfaceHitView<Real_v> *hitInfo)
  {
    // using namespace PolyconeTypes;
    if constexpr (RecordHitInfo) hitInfo->Clear();
    distance            = kInfLength;
    Vector3D<Real_v> pn = point;

    const Real_v rho2 = point.Perp2();
    int index         = polycone.GetSectionIndex(point.z());
    if (index == -1) {
      const bool onBottomCap = Abs(point.z() - Real_v(polycone.GetZAtPlane(0))) < Real_v(kTolerance);
      if (!onBottomCap || !SectionEndOwnsPoint(polycone.GetSection(0), point, rho2, false)) {
        distance = -1;
        return;
      }
      index = 0;
    } else if (index == -2) {
      const int lastIndex = polycone.GetNSections() - 1;
      const bool onTopCap = Abs(point.z() - Real_v(polycone.GetZAtPlane(polycone.GetNSections()))) < Real_v(kTolerance);
      if (!onTopCap || !SectionEndOwnsPoint(polycone.GetSection(lastIndex), point, rho2, true)) {
        distance = -1;
        return;
      }
      index = lastIndex;
    }

    // Only points close to a shared z plane need cross-section ownership
    // resolution. Ordinary starts are validated by the section DistanceToOut.
    if (index + 1 < polycone.GetNSections() &&
        Abs(point.z() - Real_v(polycone.GetZAtPlane(index + 1))) < Real_v(kTolerance)) {
      if (!ResolveSharedZSection(polycone, point, dir, rho2, index, index + 1, index)) {
        distance = -1;
        return;
      }
    } else if (index > 0 && Abs(point.z() - Real_v(polycone.GetZAtPlane(index))) < Real_v(kTolerance)) {
      if (!ResolveSharedZSection(polycone, point, dir, rho2, index - 1, index, index)) {
        distance = -1;
        return;
      }
    }

    const bool phiLimited = polycone.fDeltaPhi < kTwoPi;
    Real_v phiDistance(kInfLength);
    [[maybe_unused]] bool phiStartHit = false;
    if (phiLimited) {
      bool unused = false, outsidePhi = false;
      polycone.fPhiWedge.template GenericKernelForContainsAndInside<Real_v, true>(point, unused, outsidePhi);
      if (outsidePhi) {
        distance = -1;
        return;
      }
      bool leavingPhi = false;
      if (Abs(polycone.fPhiWedge.SafetyToOut(point)) < Real_v(kTolerance)) {
        if constexpr (RecordHitInfo) {
          const bool leavingPhiStart =
              polycone.fPhiWedge.template IsPointOnSurfaceAndMovingOut<Real_v, true, true>(point, dir);
          const bool leavingPhiEnd =
              polycone.fPhiWedge.template IsPointOnSurfaceAndMovingOut<Real_v, false, true>(point, dir);
          leavingPhi = leavingPhiStart || leavingPhiEnd;
          if (leavingPhi) phiStartHit = leavingPhiStart;
        } else {
          leavingPhi = polycone.fPhiWedge.template IsPointOnSurfaceAndMovingOut<Real_v, true, true>(point, dir) ||
                       polycone.fPhiWedge.template IsPointOnSurfaceAndMovingOut<Real_v, false, true>(point, dir);
        }
      }
      if (leavingPhi) {
        phiDistance = Real_v(0.);
      } else {
        Real_v phiDistance1(kInfLength), phiDistance2(kInfLength);
        polycone.fPhiWedge.DistanceToOut(point, dir, phiDistance1, phiDistance2);
        if (phiDistance1 < kInfLength &&
            !polycone.fPhiWedge.template IsOnSurfaceGeneric<Real_v, true>(point + phiDistance1 * dir))
          phiDistance1 = kInfLength;
        if (phiDistance2 < kInfLength &&
            !polycone.fPhiWedge.template IsOnSurfaceGeneric<Real_v, false>(point + phiDistance2 * dir))
          phiDistance2 = kInfLength;
        phiDistance = Min(phiDistance1, phiDistance2);
        if constexpr (RecordHitInfo) phiStartHit = phiDistance1 <= phiDistance2;
      }
    }

    Precision totalDistance = 0.;
    Precision dist;
    int increment = (dir.z() > 0) ? 1 : -1;
    if (Abs(dir.z()) < Real_v(kTolerance)) increment = 0;

    int istep                              = 0;
    bool wrongSideStart                    = false;
    [[maybe_unused]] bool selectedPhi      = false;
    [[maybe_unused]] bool finalHitZ        = false;
    [[maybe_unused]] bool globalEndCap     = false;
    [[maybe_unused]] int finalSectionIndex = index;
    do {
      const PolyconeSection &section = polycone.GetSection(index);

      pn = point + totalDistance * dir;
      pn.z() -= section.fShift;
      istep++;

      bool hitZSurface = false;
      // Section hand-off needs only the existing z-hit protocol. Surface-hit
      // classification is deferred until the final section has been selected.
      if (increment != 0) {
        hitZSurface = SectionImpl_t::template DistanceToOutWithZHit<Real_v>(section.fSolid, pn, dir, stepMax, dist);
      } else {
        SectionImpl_t::template DistanceToOut<Real_v>(section.fSolid, pn, dir, stepMax, dist);
      }
      if (dist == -1) {
        wrongSideStart = (istep == 1);
        break;
      }

      if (phiDistance <= Real_v(totalDistance) + dist + Real_v(kTolerance)) {
        totalDistance = phiDistance;
        if constexpr (RecordHitInfo) {
          selectedPhi       = true;
          finalSectionIndex = index;
        }
        break;
      }

      if (increment != 0 && dist < kInfLength && !hitZSurface) {
        totalDistance += dist;
        if constexpr (RecordHitInfo) finalSectionIndex = index;
        break;
      }

      totalDistance += dist;
      if (increment == 0) {
        if constexpr (RecordHitInfo) finalSectionIndex = index;
        break;
      }

      const int nextIndex = index + increment;
      if (nextIndex < 0 || nextIndex >= polycone.GetNSections()) {
        if constexpr (RecordHitInfo) {
          finalSectionIndex = index;
          finalHitZ         = hitZSurface;
          globalEndCap      = hitZSurface;
        }
        break;
      }
      if (!hitZSurface) {
        if constexpr (RecordHitInfo) finalSectionIndex = index;
        break;
      }

      const Vector3D<Real_v> hitPoint = point + totalDistance * dir;
      const bool nextOwnsHit =
          SectionEndOwnsForwardInterval(polycone.GetSection(nextIndex), hitPoint, hitPoint.Perp2(), increment < 0, dir);
      if (!nextOwnsHit) {
        if constexpr (RecordHitInfo) {
          finalSectionIndex = index;
          finalHitZ         = true;
        }
        break;
      }

      index += increment;
    } while (increment != 0 && index >= 0 && index < polycone.GetNSections());

    distance = wrongSideStart ? Real_v(-1.) : totalDistance;

    if constexpr (RecordHitInfo) {
      // A failed section solve has no boundary to associate with the hit view.
      if (wrongSideStart || distance == Real_v(kInfLength)) return;
      const Vector3D<Real_v> hitPoint     = point + distance * dir;
      const PolyconeSection &finalSection = polycone.GetSection(finalSectionIndex);
      if (selectedPhi) {
        if (polycone.fDeltaPhi <= kPi) {
          SetPhiExitHit(polycone, phiStartHit, hitInfo);
        } else {
          const bool lowerCap = Abs(hitPoint.z() - Real_v(polycone.GetZAtPlane(0))) < Real_v(kTolerance);
          const bool upperCap =
              Abs(hitPoint.z() - Real_v(polycone.GetZAtPlane(polycone.GetNSections()))) < Real_v(kTolerance);
          if (lowerCap || upperCap)
            SectionImpl_t::SetZExitHit(hitInfo, upperCap);
          else if (finalSection.fConvex && IsIncidentOnOuterSurface(finalSection, hitPoint))
            SetOuterExitHit(finalSection, hitPoint, hitInfo);
          else
            SetPhiExitHit(polycone, phiStartHit, hitInfo);
        }
      } else if (finalHitZ) {
        if (globalEndCap)
          SectionImpl_t::SetZExitHit(hitInfo, increment > 0);
        else if (finalSection.fConvex && IsIncidentOnOuterSurface(finalSection, hitPoint))
          SetOuterExitHit(finalSection, hitPoint, hitInfo);
      } else if (finalSection.fConvex && IsIncidentOnOuterSurface(finalSection, hitPoint)) {
        SetOuterExitHit(finalSection, hitPoint, hitInfo);
      }
    }

    return;
  }

  /**
   * @brief Compute the inside-to-exit distance and optionally record a supporting exit normal.
   * @param polycone Polycone geometry and section-support cache.
   * @param point Starting point in the polycone frame.
   * @param dir Normalized ray direction.
   * @param stepMax Caller-provided limit; existing distance behavior is unchanged.
   * @param[out] distance Exit distance, or `-1` for a rejected outside point.
   * @param hitInfo Optional hit view; a non-null pointer requests `leftBehind` and, when supplied by the view,
   *                the corresponding unit normal.
   * @details Non-supporting exits leave normal storage untouched. Section-local
   * z crossings remain internal handoffs, not supporting global endcaps.
   */
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void DistanceToOut(UnplacedStruct_t const &polycone,
                                                                         Vector3D<Real_v> const &point,
                                                                         Vector3D<Real_v> const &dir,
                                                                         Real_v const &stepMax, Real_v &distance,
                                                                         SurfaceHitView<Real_v> *hitInfo = nullptr)
  {
    if (hitInfo)
      DistanceToOutImpl<Real_v, true>(polycone, point, dir, stepMax, distance, hitInfo);
    else
      DistanceToOutImpl<Real_v, false>(polycone, point, dir, stepMax, distance, nullptr);
  }

  /**
   * @brief Compute the outside safety for the polycone.
   *
   * The method starts from the seeded section-cone safety and refines it with a
   * local search on the relevant `(r,z)` contour. This is required for sharp
   * jumps and repeated-`z` connectors where the nearest polycone boundary is not
   * owned by the seeded section cone.
   */
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SafetyToIn(UnplacedStruct_t const &polycone,
                                                                      Vector3D<Real_v> const &point, Real_v &safety)
  {

    using namespace ConeUtilities;

    Vector3D<Real_v> p = point;
    int index          = polycone.GetSectionIndex(p.z());

    bool needZ = false;
    if (index < 0) {
      needZ = true;
      if (index == -1) index = 0;
      if (index == -2) index = polycone.GetNSections() - 1;
    }

    PolyconeSection const &sec = polycone.GetSection(index);
    Vector3D<Real_v> localp    = p - Vector3D<Precision>(0, 0, sec.fShift);
    SectionImpl_t::template SafetyToIn<Real_v>(sec.fSolid, localp, safety);

    // Match the conservative cone-safety ordering: phi is needed to classify
    // a point that is radially inside but outside the wedge. When the radial
    // union already supplies a non-negative lower bound, evaluating the wedge
    // cannot make that bound safer and was redundant in the legacy kernel too.
    if (polycone.fDeltaPhi < kTwoPi && safety < Real_v(0.)) safety = Max(safety, polycone.fPhiWedge.SafetyToIn(point));

    if (safety < kTolerance) return;

    Real_v best = safety;

    const Real_v rho2 = p.x() * p.x() + p.y() * p.y();

    Real_v clampedLocalZ     = localp.z();
    clampedLocalZ            = Max(Real_v(-sec.fSolid.fDz), Min(Real_v(sec.fSolid.fDz), clampedLocalZ));
    const Real_v outerRadius = GetRadiusOfConeAtPoint<Real_v, false>(sec.fSolid, clampedLocalZ);
    Real_v innerRadius(0.);
    if (sec.fSolid.fRmin1 > 0. || sec.fSolid.fRmin2 > 0.) {
      innerRadius = GetRadiusOfConeAtPoint<Real_v, true>(sec.fSolid, clampedLocalZ);
    }

    const Vector<Vector3D<Precision>> *contour = nullptr;
    const Real_v outerLimit                    = outerRadius + Real_v(kTolerance);
    if (rho2 > outerLimit * outerLimit) {
      contour = &polycone.fRMaxTwoDVec;
    } else if (sec.fSolid.fRmin1 > 0. || sec.fSolid.fRmin2 > 0.) {
      const Real_v innerLimit = innerRadius - Real_v(kTolerance);
      if (innerLimit > Real_v(0.) && rho2 < innerLimit * innerLimit) {
        contour = &polycone.fRMinTwoDVec;
      }
    } else {
      return;
    }
    if (!contour) return;

    const int npts = static_cast<int>(contour->size());
    if (npts < 2) return;

    // For outside safety the nearest polycone boundary can be a profile jump or
    // repeated-z connector that the seeded cone section does not own. Search
    // the relevant rz contour locally, pruning by the monotone dz^2 bound.
    Real_v bestSq    = (best < Real_v(kInfLength)) ? best * best : Real_v(kInfLength);
    const Real_v rho = Sqrt(rho2);

    int seedSegment = 2 * index;
    if (needZ && p.z() < polycone.fZs[0]) seedSegment = 0;
    if (needZ && p.z() > polycone.fZs[polycone.fZs.size() - 1]) seedSegment = npts - 2;
    seedSegment = Max(0, Min(npts - 2, seedSegment));

    for (int i = seedSegment; i >= 0; --i) {
      const auto &a        = (*contour)[i];
      const auto &b        = (*contour)[i + 1];
      const Real_v minDzSq = SegmentMinDeltaZSquared(p.z(), a, b);
      if (bestSq < Real_v(kInfLength) && minDzSq >= bestSq) break;

      const Real_v distSq = DistanceSquaredToRZSegment(rho, p.z(), a, b);
      if (distSq < bestSq) bestSq = distSq;
    }

    for (int i = seedSegment + 1; i < (npts - 1); ++i) {
      const auto &a        = (*contour)[i];
      const auto &b        = (*contour)[i + 1];
      const Real_v minDzSq = SegmentMinDeltaZSquared(p.z(), a, b);
      if (bestSq < Real_v(kInfLength) && minDzSq >= bestSq) break;

      const Real_v distSq = DistanceSquaredToRZSegment(rho, p.z(), a, b);
      if (distSq < bestSq) bestSq = distSq;
    }

    const Real_v contourSafety = Sqrt(bestSq);
    if (contourSafety < safety) safety = contourSafety;
    return;
  }

  /**
   * @brief Compute the inside safety for the polycone.
   *
   * The nearest exit is found from the closed `(r,z)` contour and then
   * minimized with the phi-wedge safety when the solid is phi-limited.
   */
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SafetyToOut(UnplacedStruct_t const &polycone,
                                                                       Vector3D<Real_v> const &point, Real_v &safety)
  {
    const Real_v z = point.z();
    safety = Min(z - Real_v(polycone.GetZAtPlane(0)), Real_v(polycone.GetZAtPlane(polycone.GetNSections())) - z);
    // The terminal z bounds already give a negative safety outside their tolerance band.
    if (safety < -Real_v(kTolerance)) return;

    int index   = polycone.GetSectionIndex(z);
    bool compIn = false, compOut = false;
    // Retain radial rejection for the unsigned contour distance, reusing the local section.
    if (index >= 0 && Abs(z - Real_v(polycone.GetZAtPlane(index))) > Real_v(kTolerance) &&
        Abs(z - Real_v(polycone.GetZAtPlane(index + 1))) > Real_v(kTolerance))
      GenericKernelForASection<Real_v, true>(polycone, index, point, compIn, compOut);
    else
      GenericKernelForContainsAndInside<Real_v, true>(polycone, point, compIn, compOut);
    if (compOut) {
      safety = Real_v(-1.);
      return;
    }

    const Real_v rho = point.Perp();
    Real_v bestSq    = Real_v(kInfLength);
    const auto visit = [&](Vector3D<Precision> const &a, Vector3D<Precision> const &b) {
      if ((a.x() == 0 && b.x() == 0) || (a.x() == b.x() && a.y() == b.y())) return;
      const Real_v distanceSq = DistanceSquaredToRZSegment(rho, z, a, b);
      if (distanceSq < bestSq) bestSq = distanceSq;
    };

    const int npts = polycone.fRMaxTwoDVec.size();
    if (npts == 2) {
      // A single section has no contour search to prune.
      for (unsigned int i = 0; i < polycone.fTwoDVec.size(); ++i)
        visit(polycone.fTwoDVec[i], polycone.fTwoDVec[i + 1 == polycone.fTwoDVec.size() ? 0 : i + 1]);
    } else {
      // Only terminal caps are physical boundaries; intermediate section caps are not exits.
      // Preserve the closed-contour endpoint order to avoid roundoff at exact cap surface points.
      visit(polycone.fRMinTwoDVec[0], polycone.fRMaxTwoDVec[0]);
      visit(polycone.fRMaxTwoDVec[npts - 1], polycone.fRMinTwoDVec[npts - 1]);
      if (index < 0) index = index == -1 ? 0 : polycone.GetNSections() - 1;
      const int seedSegment = 2 * index;

      // Both profiles have monotone z, including radial jumps and repeated-z connectors.
      // Search outwards from the local section until z separation exceeds the best distance.
      for (int profile = 0; profile < 2; ++profile) {
        const auto &contour = profile == 0 ? polycone.fRMaxTwoDVec : polycone.fRMinTwoDVec;
        for (int i = seedSegment; i >= 0; --i) {
          if (SegmentMinDeltaZSquared(z, contour[i], contour[i + 1]) > bestSq) break;
          visit(contour[i], contour[i + 1]);
        }
        for (int i = seedSegment + 1; i < npts - 1; ++i) {
          if (SegmentMinDeltaZSquared(z, contour[i], contour[i + 1]) > bestSq) break;
          visit(contour[i], contour[i + 1]);
        }
      }
    }

    safety = Sqrt(bestSq);

    if (polycone.fDeltaPhi < 2 * kPi) {
      Real_v safetyPhi = polycone.fPhiWedge.SafetyToOut(point);
      if (safetyPhi < safety) safety = safetyPhi;
    }
  }
};
} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif // VECGEOM_VOLUMES_KERNEL_polyconeIMPLEMENTATION_H_
