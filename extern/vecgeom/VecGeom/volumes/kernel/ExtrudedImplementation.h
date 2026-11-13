//===-- kernel/ExtrudedImplementation.h ----------------------------------*- C++ -*-===//
//===--------------------------------------------------------------------------===//
/// @file ExtrudedImplementation.h
/// @brief Navigation dispatcher for simple and tessellated extruded solids.
/// @author mihaela.gheata@cern.ch

#ifndef VECGEOM_VOLUMES_KERNEL_EXTRUDEDIMPLEMENTATION_H_
#define VECGEOM_VOLUMES_KERNEL_EXTRUDEDIMPLEMENTATION_H_

#include <cstdio>
#include <VecCore/VecCore>
#include "VecGeom/base/Config.h"
#include "VecGeom/volumes/kernel/GenericKernels.h"
#include "VecGeom/base/Vector3D.h"
#include "VecGeom/management/ABBoxImplementation.h"
#include "VecGeom/volumes/SurfaceHitView.h"

#include "TessellatedImplementation.h"
#include "SExtruImplementation.h"

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE(struct ExtrudedImplementation;);
VECGEOM_DEVICE_DECLARE_CONV(struct, ExtrudedImplementation);

inline namespace VECGEOM_IMPL_NAMESPACE {

class PlacedExtruded;
class ExtrudedStruct;
class UnplacedExtruded;

/// @brief Dispatches extruded-solid navigation to SExtru or tessellated helpers.
/// @details Single-section simple extrusions use `SExtruImplementation`.
/// General multi-section extrusions use tessellated helpers, with an optional
/// per-section path when `fUseTslSections` is enabled.
struct ExtrudedImplementation {

  using PlacedShape_t    = PlacedExtruded;
  using UnplacedStruct_t = ExtrudedStruct;
  using UnplacedVolume_t = UnplacedExtruded;

private:
  /// @brief Publish a section-facet normal using its construction-time support classification.
  template <typename Real_v, typename Facet_t>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SetSectionFacetExitHit(Facet_t const &facet,
                                                                                  SurfaceHitView<Real_v> *hitInfo)
  {
    bool const leftBehind = facet.fConvex;
    hitInfo->SetLeftBehind(leftBehind);
    if (leftBehind && hitInfo->WantsNormal()) hitInfo->SetNormal(Vector3D<Real_v>(facet.fNormal));
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static int FindRaySection(UnplacedStruct_t const &extruded,
                                                                         Real_v const &pointZ, Real_v const &directionZ)
  {
    const int nseg = static_cast<int>(extruded.GetNSegments());
    int zIndex     = extruded.FindZSegment(static_cast<Precision>(pointZ));
    if (zIndex > 0 && zIndex < nseg && directionZ < Real_v(0.) &&
        Abs(pointZ - Real_v(extruded.fZPlanes[zIndex])) <= Real_v(kTolerance)) {
      --zIndex;
    }
    if (zIndex < 0) return 0;
    if (zIndex >= nseg) return nseg - 1;
    return zIndex;
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void DistanceToInSections(UnplacedStruct_t const &extruded,
                                                                                Vector3D<Real_v> const &point,
                                                                                Vector3D<Real_v> const &direction,
                                                                                Real_v const &stepMax, Real_v &distance)
  {
    const int nseg = static_cast<int>(extruded.GetNSegments());
    distance       = InfinityLength<Real_v>();
    if (nseg <= 0) return;

    const Real_v invdirz = Real_v(1.) / NonZero(direction.z());
    int zIndex           = FindRaySection(extruded, point.z(), direction.z());
    const int step       = direction.z() >= Real_v(0.) ? 1 : -1;

    for (; zIndex >= 0 && zIndex < nseg; zIndex += step) {
      // Each section owns one z interval; the first accepted interval hit is
      // the first material entry along the ray.
      const Real_v candidate =
          extruded.fTslSections[zIndex]->template DistanceToIn<false>(point, direction, invdirz, stepMax);
      if (candidate < Real_v(-kTolerance)) {
        // A tolerated surface offset divided by a shallow inward projection
        // can place the section interval just behind the ray origin. Classify
        // only this ambiguous path so ordinary outside queries stay unchanged.
        Inside_t inside;
        InsideSectionsWithRuntimeFallback<Real_v>(extruded, point, inside);
        if (inside != kOutside) {
          distance = inside == kSurface ? Real_v(0.) : -Real_v(1.);
          return;
        }
        continue;
      }
      if (candidate >= Real_v(-kTolerance) && candidate < stepMax && candidate < InfinityLength<Real_v>()) {
        distance = Max(candidate, Real_v(0.));
        return;
      }
      if (direction.z() == Real_v(0.)) return;
    }
  }

  template <typename Real_v, bool RecordHitInfo>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void DistanceToOutSections(
      UnplacedStruct_t const &extruded, Vector3D<Real_v> const &point, Vector3D<Real_v> const &direction,
      Real_v const & /*stepMax*/, Real_v &distance, SurfaceHitView<Real_v> *hitInfo)
  {
    const int nseg = static_cast<int>(extruded.GetNSegments());
    distance       = -Real_v(1.);
    if (nseg <= 0) return;
    if (point.z() < Real_v(extruded.fZPlanes[0]) - Real_v(kTolerance) ||
        point.z() > Real_v(extruded.fZPlanes[nseg]) + Real_v(kTolerance)) {
      return;
    }
    // Starts on a global z cap moving out of, or tangent to, the solid leave at zero.
    if ((Abs(point.z() - Real_v(extruded.fZPlanes[nseg])) <= Real_v(kTolerance) && direction.z() >= Real_v(0.)) ||
        (Abs(point.z() - Real_v(extruded.fZPlanes[0])) <= Real_v(kTolerance) && direction.z() <= Real_v(0.))) {
      distance = Real_v(0.);
      if constexpr (RecordHitInfo) {
        // Tangential rays have no z sign; select the cap from the start point.
        bool const upper = Abs(point.z() - Real_v(extruded.fZPlanes[nseg])) <= Real_v(kTolerance);
        hitInfo->SetLeftBehind(true);
        if (hitInfo->WantsNormal()) hitInfo->SetNormal(Vector3D<Real_v>(0., 0., upper ? 1. : -1.));
      }
      return;
    }

    const Real_v invdirz = Real_v(1.) / NonZero(direction.z());
    int zIndex           = FindRaySection(extruded, point.z(), direction.z());
    const int step       = direction.z() >= Real_v(0.) ? 1 : -1;

    for (; zIndex >= 0 && zIndex < nseg; zIndex += step) {
      // Internal z planes connect adjacent sections; only side exits from the
      // current segment or the final global cap can leave the full extruded.
      int hitFacet = -1;
      Real_v sideDistance;
      if constexpr (RecordHitInfo)
        sideDistance =
            extruded.fTslSections[zIndex]->template DistanceToOutRange<true>(point, direction, invdirz, &hitFacet);
      else
        sideDistance = extruded.fTslSections[zIndex]->DistanceToOutRange(point, direction, invdirz);
      if (sideDistance < InfinityLength<Real_v>()) {
        // The start is required to be inside or on surface (checked by the
        // public wrong-side guard when enabled). A tolerated spatial offset
        // can yield a larger backward ray parameter at shallow incidence.
        distance = Max(sideDistance, Real_v(0.));
        if constexpr (RecordHitInfo) {
          if (hitFacet >= 0) {
            auto const &facet = extruded.fTslSections[zIndex]->GetFacet(hitFacet);
            SetSectionFacetExitHit<Real_v>(facet, hitInfo);
          }
        }
        return;
      }
      if (direction.z() == Real_v(0.)) {
        distance = InfinityLength<Real_v>();
        return;
      }
      const int nextIndex = zIndex + step;
      if (nextIndex < 0 || nextIndex >= nseg) {
        const Real_v endZ = step > 0 ? Real_v(extruded.fZPlanes[nseg]) : Real_v(extruded.fZPlanes[0]);
        distance          = Max((endZ - point.z()) / direction.z(), Real_v(0.));
        if constexpr (RecordHitInfo) {
          hitInfo->SetLeftBehind(true);
          if (hitInfo->WantsNormal()) hitInfo->SetNormal(Vector3D<Real_v>(0., 0., step > 0 ? 1. : -1.));
        }
        return;
      }
    }
  }

  /// @brief Dispatch an exit query while compiling hit recording completely out of the ordinary path.
  template <typename Real_v, bool RecordHitInfo>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void DistanceToOutImpl(UnplacedStruct_t const &extruded,
                                                                             Vector3D<Real_v> const &point,
                                                                             Vector3D<Real_v> const &direction,
                                                                             Real_v const &stepMax, Real_v &distance,
                                                                             SurfaceHitView<Real_v> *hitInfo)
  {
    if constexpr (RecordHitInfo) hitInfo->Clear();

    if (extruded.fIsSxtru) {
      SExtruImplementation::DistanceToOut<Real_v>(extruded.fSxtruHelper, point, direction, stepMax, distance,
                                                  RecordHitInfo ? hitInfo : nullptr);
      return;
    }

    if (extruded.fUseTslSections) {
#ifndef VECGEOM_NO_WRONGSIDE_CONV
      Inside_t inside;
      InsideSectionsWithRuntimeFallback<Real_v>(extruded, point, inside);
      // Section helpers are accelerator kernels; public wrong-side conventions
      // must be decided before asking them for geometry exits.
      if (inside == kOutside) {
        distance = Real_v(-1.);
        return;
      }
#endif
      DistanceToOutSections<Real_v, RecordHitInfo>(extruded, point, direction, stepMax, distance, hitInfo);
      return;
    }

    TessellatedImplementation::DistanceToOut<Real_v>(extruded.fTslRuntimeHelper, point, direction, stepMax, distance,
                                                     RecordHitInfo ? hitInfo : nullptr);
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SafetyToInSections(UnplacedStruct_t const &extruded,
                                                                              Vector3D<Real_v> const &point,
                                                                              Real_v &safety)
  {
    const int nseg = static_cast<int>(extruded.GetNSegments());
    safety         = InfinityLength<Real_v>();
    for (int iseg = 0; iseg < nseg; ++iseg)
      safety = Min(safety, Real_v(extruded.fTslSections[iseg]->SafetyToIn(point)));
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SafetyToOutSections(UnplacedStruct_t const &extruded,
                                                                               Vector3D<Real_v> const &point,
                                                                               Real_v &safety)
  {
    const int nseg = static_cast<int>(extruded.GetNSegments());
    safety         = -Real_v(1.);
    if (nseg <= 0) return;
    if (point.z() < Real_v(extruded.fZPlanes[0]) - Real_v(kTolerance) ||
        point.z() > Real_v(extruded.fZPlanes[nseg]) + Real_v(kTolerance)) {
      return;
    }

    const int zIndex     = FindRaySection(extruded, point.z(), Real_v(1.));
    const Real_v zCapSaf = Min(point.z() - Real_v(extruded.fZPlanes[0]), Real_v(extruded.fZPlanes[nseg]) - point.z());
    safety               = Min(extruded.fTslSections[zIndex]->LateralSafetyToOut(point), zCapSaf);
    // The closest side point can lie on a kink shared with an adjacent section;
    // internal z planes are still not exits, so only lateral neighbor planes
    // are allowed to reduce the conservative safety.
    if (zIndex > 0) {
      const Real_v neighborSafety = extruded.fTslSections[zIndex - 1]->LateralSafetyToOut(point);
      if (neighborSafety > Real_v(0.)) safety = Min(safety, neighborSafety);
    }
    if (zIndex + 1 < nseg) {
      const Real_v neighborSafety = extruded.fTslSections[zIndex + 1]->LateralSafetyToOut(point);
      if (neighborSafety > Real_v(0.)) safety = Min(safety, neighborSafety);
    }
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void InsideSectionsWithRuntimeFallback(
      UnplacedStruct_t const &extruded, Vector3D<Real_v> const &point, Inside_t &inside)
  {
    Inside<Real_v>(extruded, point, inside);
    if (inside != kOutside) return;

    // Section classification is the cheap first pass, but wrong-side sentinels
    // must not reject tolerated tessellated surface points owned by the runtime
    // helper.
    Inside_t runtimeInside;
    TessellatedImplementation::Inside<Real_v>(extruded.fTslRuntimeHelper, point, runtimeInside);
    if (runtimeInside != kOutside) inside = runtimeInside;
  }

  /// @brief Share polygon projection while retaining each query's boundary convention.
  template <typename Real_v, bool ForInside>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void ClassifyRuntime(UnplacedStruct_t const &extruded,
                                                                           Vector3D<Real_v> const &point,
                                                                           Inside_t &inside)
  {
    inside                  = kOutside;
    const auto &runtime     = extruded.fTslRuntimeHelper;
    constexpr Real_v margin = ForInside ? Real_v(kHalfTolerance) : Real_v(0.);
    bool contains           = false;
    ABBoxImplementation::ABBoxContainsKernel(runtime.fMinExtent - Vector3D<Real_v>(margin),
                                             runtime.fMaxExtent + Vector3D<Real_v>(margin), point, contains);
    if (!contains) return;

    // A Contains cutoff displacement hidden by coordinate roundoff cannot
    // resolve a surface point; retain the original classifier without probing.
    constexpr Real_v rayCutoff = Real_v(1.e-8);
    if (extruded.fProjectionSafetySqr > 0. && (ForInside || extruded.fProjectionZTolerance < rayCutoff)) {
      // Contains counts crossings after 1.e-8 along its unit test ray. Classify
      // beyond that cutoff, retaining the original ray test when the displaced
      // point cannot be separated from a boundary by the cached uncertainty band.
      const Vector3D<Real_v> query = ForInside ? point : point + rayCutoff * runtime.fTestDir;
      if constexpr (!ForInside) {
        if (query.z() < extruded.fSections[0].fOrigin.z() - extruded.fProjectionZTolerance ||
            query.z() > extruded.fSections[1].fOrigin.z() + extruded.fProjectionZTolerance)
          return;
      }
      if (query.z() - extruded.fSections[0].fOrigin.z() > extruded.fProjectionZTolerance &&
          extruded.fSections[1].fOrigin.z() - query.z() > extruded.fProjectionZTolerance) {
        const auto &lower     = extruded.fSections[0];
        const auto &upper     = extruded.fSections[1];
        const Real_v fraction = (query.z() - lower.fOrigin.z()) / (upper.fOrigin.z() - lower.fOrigin.z());
        const Real_v scale    = (1. - fraction) * lower.fScale + fraction * upper.fScale;
        const Vector3D<Real_v> projected(
            (query.x() - ((1. - fraction) * lower.fOrigin.x() + fraction * upper.fOrigin.x())) / scale,
            (query.y() - ((1. - fraction) * lower.fOrigin.y() + fraction * upper.fOrigin.y())) / scale, 0.);
        int closest;
        // A cross-section edge moves at most V per unit z. Its 2D separation d
        // bounds 3D lateral separation below by d/sqrt(1+V^2). The cached band
        // includes that slope bound, both endpoint scales and projection roundoff.
        if (extruded.fPolygon.SafetySqr(projected, closest) > extruded.fProjectionSafetySqr) {
          inside = extruded.fPolygon.Contains(projected) ? kInside : kOutside;
          return;
        }
        if constexpr (ForInside) {
          if (closest >= 0) {
            // Polygon edge i joins vertices i and i-1. Construction stores the lower
            // cap's N-2 triangles first, followed by two lateral triangles per edge.
            const size_t nvertices = extruded.fPolygon.GetNVertices();
            const size_t side      = closest == 0 ? nvertices - 1 : size_t(closest - 1);
            const size_t facet     = nvertices - 2 + 2 * side;
            if (TessellatedImplementation::IsOnFiniteFacet(runtime.fFacets[facet], point) ||
                TessellatedImplementation::IsOnFiniteFacet(runtime.fFacets[facet + 1], point)) {
              inside = kSurface;
              return;
            }
          }
        }
      }
    }
    TessellatedImplementation::ClassifyInside<Real_v, ForInside>(runtime, point, inside);
  }

public:
  /// @brief Test whether a local point is inside or on the tolerated surface.
  /// @details Supported concave extrusions project the point beyond the existing
  /// tessellated ray cutoff. Unresolved boundaries retain the original crossing test.
  /// @tparam Real_v Floating-point scalar type.
  /// @param extruded Runtime extruded data.
  /// @param point Local point to classify.
  /// @param[out] inside True when @p point is contained by the selected helper.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void Contains(UnplacedStruct_t const &extruded,
                                                                    Vector3D<Real_v> const &point, bool &inside)
  {
    inside = false;
    if (extruded.fIsSxtru) {
      SExtruImplementation::Contains<Real_v>(extruded.fSxtruHelper, point, inside);
      return;
    }

    if (extruded.fUseTslSections) {
      // Find the Z section
      const int nseg = (int)extruded.GetNSegments();
      if (point.z() < Real_v(extruded.fZPlanes[0]) - Real_v(kTolerance) ||
          point.z() > Real_v(extruded.fZPlanes[nseg]) + Real_v(kTolerance)) {
        return;
      }
      int zIndex = extruded.FindZSegment(point[2]);
      if ((zIndex < 0) || (zIndex > nseg)) return;
      inside = extruded.fTslSections[Min(zIndex, nseg - 1)]->Contains(point);
      return;
    }
    Inside_t classification;
    ClassifyRuntime<Real_v, false>(extruded, point, classification);
    inside = classification == kInside;
  }

  /// @brief Classify a local point as inside, outside, or surface.
  /// @details Dispatches to the simple-extruded helper, the optional
  /// per-section tessellated helper, or the runtime tessellated helper. The
  /// per-section path promotes points on the first or last z section to
  /// `kSurface` when they are within tolerance of the section plane.
  /// Two-section concave extrusions use polygon projection when certified away
  /// from every face; the runtime helper retains all near-surface decisions.
  /// @tparam Real_v Floating-point scalar type.
  /// @param extruded Runtime extruded data.
  /// @param point Local point to classify.
  /// @param[out] inside Set to `kInside`, `kOutside`, or `kSurface`.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void Inside(UnplacedStruct_t const &extruded,
                                                                  Vector3D<Real_v> const &point, Inside_t &inside)
  {
    inside = EInside::kOutside;

    if (extruded.fIsSxtru) {
      SExtruImplementation::Inside<Real_v>(extruded.fSxtruHelper, point, inside);
      return;
    }

    if (extruded.fUseTslSections) {
      const int nseg = (int)extruded.GetNSegments();
      if (point.z() < Real_v(extruded.fZPlanes[0]) - Real_v(kTolerance) ||
          point.z() > Real_v(extruded.fZPlanes[nseg]) + Real_v(kTolerance)) {
        return;
      }
      int zIndex = extruded.FindZSegment(point[2]);
      if ((zIndex < 0) || (zIndex > nseg)) return;
      inside = extruded.fTslSections[Min(zIndex, nseg - 1)]->Inside(point);
      if (inside == EInside::kOutside) return;
      if (inside == EInside::kInside) {
        // Need to check if point on Z section
        if (((zIndex == 0) || (zIndex == nseg)) && Abs(point[2] - extruded.fZPlanes[zIndex]) < kTolerance) {
          inside = EInside::kSurface;
        }
      } else {
        inside = EInside::kSurface;
      }
      return;
    }
    ClassifyRuntime<Real_v, true>(extruded, point, inside);
  }

  /// @brief Compute distance from an outside or surface point to enter the solid.
  /// @details Simple extrusions are delegated to `SExtruImplementation`.
  /// General extrusions use the tessellated runtime helper unless the optional
  /// section accelerator is compiled in.
  /// @tparam Real_v Floating-point scalar type.
  /// @param extruded Runtime extruded data.
  /// @param point Local start point.
  /// @param direction Unit local direction.
  /// @param stepMax Maximum distance to consider.
  /// @param[out] distance Entry distance, `-1` for wrong-side starts in public
  /// helper paths, or `kInfLength` when no entry is found.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void DistanceToIn(UnplacedStruct_t const &extruded,
                                                                        Vector3D<Real_v> const &point,
                                                                        Vector3D<Real_v> const &direction,
                                                                        Real_v const &stepMax, Real_v &distance)
  {
    if (extruded.fIsSxtru)
      SExtruImplementation::DistanceToIn<Real_v>(extruded.fSxtruHelper, point, direction, stepMax, distance);
    else if (extruded.fUseTslSections) {
      DistanceToInSections<Real_v>(extruded, point, direction, stepMax, distance);
      if (distance == InfinityLength<Real_v>()) {
        TessellatedImplementation::DistanceToIn<Real_v>(extruded.fTslRuntimeHelper, point, direction, stepMax,
                                                        distance);
      }
    } else
      TessellatedImplementation::DistanceToIn<Real_v>(extruded.fTslRuntimeHelper, point, direction, stepMax, distance);
  }

  /// @brief Compute distance from an inside or surface point to leave the solid.
  /// @tparam Real_v Floating-point scalar type.
  /// @param extruded Runtime extruded data.
  /// @param point Local start point.
  /// @param direction Unit local direction.
  /// @param stepMax Maximum distance to consider.
  /// @param[out] distance Exit distance, `-1` for wrong-side starts in public
  /// helper paths, or the tessellated helper miss convention for finite
  /// `stepMax` limits.
  /// @param[out] hitInfo Optional supporting exit normal; false validity leaves normal storage untouched.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void DistanceToOut(UnplacedStruct_t const &extruded,
                                                                         Vector3D<Real_v> const &point,
                                                                         Vector3D<Real_v> const &direction,
                                                                         Real_v const &stepMax, Real_v &distance,
                                                                         SurfaceHitView<Real_v> *hitInfo = nullptr)
  {
    if (hitInfo)
      DistanceToOutImpl<Real_v, true>(extruded, point, direction, stepMax, distance, hitInfo);
    else
      DistanceToOutImpl<Real_v, false>(extruded, point, direction, stepMax, distance, nullptr);
  }

  /// @brief Compute safety from an outside point to enter the solid.
  /// @tparam Real_v Floating-point scalar type.
  /// @param extruded Runtime extruded data.
  /// @param point Local point.
  /// @param[out] safety Conservative distance to the nearest entry boundary,
  /// with public helper sentinels supplied by the selected implementation.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SafetyToIn(UnplacedStruct_t const &extruded,
                                                                      Vector3D<Real_v> const &point, Real_v &safety)
  {
    if (extruded.fIsSxtru)
      SExtruImplementation::SafetyToIn<Real_v>(extruded.fSxtruHelper, point, safety);
    else if (extruded.fUseTslSections) {
#ifndef VECGEOM_NO_WRONGSIDE_CONV
      Inside_t inside;
      InsideSectionsWithRuntimeFallback<Real_v>(extruded, point, inside);
      // Entry safety is wrong-sided from material and exactly zero on the
      // tolerated surface, independent of the section accelerator.
      if (inside == kInside) {
        safety = Real_v(-1.);
        return;
      }
      if (inside == kSurface) {
        safety = Real_v(0.);
        return;
      }
#endif
      SafetyToInSections<Real_v>(extruded, point, safety);
    } else
      TessellatedImplementation::SafetyToIn<Real_v>(extruded.fTslRuntimeHelper, point, safety);
  }

  /// @brief Compute safety from an inside point to leave the solid.
  /// @tparam Real_v Floating-point scalar type.
  /// @param extruded Runtime extruded data.
  /// @param point Local point.
  /// @param[out] safety Conservative distance to the nearest exit boundary,
  /// with public helper sentinels supplied by the selected implementation.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SafetyToOut(UnplacedStruct_t const &extruded,
                                                                       Vector3D<Real_v> const &point, Real_v &safety)
  {
    if (extruded.fIsSxtru)
      SExtruImplementation::SafetyToOut<Real_v>(extruded.fSxtruHelper, point, safety);
    else if (extruded.fUseTslSections) {
#ifndef VECGEOM_NO_WRONGSIDE_CONV
      Inside_t inside;
      InsideSectionsWithRuntimeFallback<Real_v>(extruded, point, inside);
      // Exit safety is only defined from material; outside starts keep the
      // public wrong-side sentinel before using section-local distances.
      if (inside == kOutside) {
        safety = Real_v(-1.);
        return;
      }
      if (inside == kSurface) {
        safety = Real_v(0.);
        return;
      }
#endif
      SafetyToOutSections<Real_v>(extruded, point, safety);
    } else
      TessellatedImplementation::SafetyToOut<Real_v>(extruded.fTslRuntimeHelper, point, safety);
  }

  /// @brief Compute a surface normal from the selected extruded helper.
  /// @tparam Real_v Floating-point scalar type.
  /// @param extruded Runtime extruded data.
  /// @param point Local surface point.
  /// @param[out] valid True when a surface normal could be assigned.
  /// @return Outward normal, or a fallback vector when @p valid is false.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static Vector3D<Real_v> NormalKernel(UnplacedStruct_t const &extruded,
                                                                                    Vector3D<Real_v> const &point,
                                                                                    bool &valid)
  {
    // Computes the normal on a surface and returns it as a unit vector
    if (extruded.fIsSxtru) return SExtruImplementation::NormalKernel<Real_v>(extruded.fSxtruHelper, point, valid);
    return TessellatedImplementation::NormalKernel<Real_v>(extruded.fTslRuntimeHelper, point, valid);
  }

}; // end ExtrudedImplementation

} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif // VECGEOM_VOLUMES_KERNEL_EXTRUDEDIMPLEMENTATION_H_
