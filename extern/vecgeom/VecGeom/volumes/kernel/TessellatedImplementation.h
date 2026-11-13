//===-- kernel/TessellatedImplementation.h ----------------------------------*- C++ -*-===//
//===--------------------------------------------------------------------------===//
/// @file TessellatedImplementation.h
/// @brief Navigation kernels for tessellated runtime solids.
/// @author mihaela.gheata@cern.ch, sandro.wenzel@cern.ch

#ifndef VECGEOM_VOLUMES_KERNEL_TESSELLATEDIMPLEMENTATION_H_
#define VECGEOM_VOLUMES_KERNEL_TESSELLATEDIMPLEMENTATION_H_

#include "VecGeom/base/Config.h"
#include "VecGeom/base/Vector3D.h"
#include "VecGeom/management/ABBoxImplementation.h"
#include "VecGeom/volumes/SurfaceHitView.h"
#include "VecGeom/volumes/TessellatedStruct.h"
#include "VecGeom/volumes/kernel/GenericKernels.h"
#include <VecCore/VecCore>

#include <cstdio>

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE(struct TessellatedImplementation;);
VECGEOM_DEVICE_DECLARE_CONV(struct, TessellatedImplementation);

inline namespace VECGEOM_IMPL_NAMESPACE {

class PlacedTessellated;
template <size_t NVERT, typename T>
class TessellatedStruct;
class UnplacedTessellated;

/// IndexedBVH leaves can contain several facets. Retest each conservatively rounded facet AABB so
/// a ray crossing the leaf bounds does not trigger expensive distance work for unrelated facets.
/// The outward-rounded boxes make this extra rejection conservative.
static constexpr bool kCheckFacetBB = true;

/// @brief Implements tessellated-solid navigation using BVH facet queries.
/// @details Point classification uses the first oriented crossing of a fixed
/// test ray, with half-open parity for ambiguous shared-edge/vertex crossings.
/// Distance and safety helpers query the runtime facet BVH and leave public
/// wrong-side sentinels to the wrapper methods.
struct TessellatedImplementation {

  using PlacedShape_t    = PlacedTessellated;
  using UnplacedStruct_t = TessellatedRuntimeStruct<Precision>;
  using UnplacedVolume_t = UnplacedTessellated;

  /**
   * @brief Traverse ray-intersected facet boxes near-first.
   * @details The visitor receives a mutable traversal bound initialized from
   * @p step. Reducing it prunes later nodes and boxes; returning true skips the
   * remaining facets in the current leaf, or ends traversal when stop_at_hit.
   * @tparam stop_at_hit End the entire traversal on a true visitor result.
   * @tparam check_leaf_bb Retest each facet AABB within an accepted leaf.
   * @tparam Real_v Query scalar type.
   * @tparam LeafFunction Facet visitor type.
   * @param bvh Facet hierarchy.
   * @param point Local ray origin.
   * @param direction Local ray direction.
   * @param step Initial ray bound.
   * @param visit Exact facet operation receiving `(facetId, traversalBound, approach)`;
   *              it may reduce the ray bound by reference. The approach is a
   *              conservative distance to the facet box (zero when not tested).
   */
  template <bool check_leaf_bb, bool stop_at_hit = false, typename Real_v, typename LeafFunction>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void VisitRayBoxes(FacetBVH_t const &bvh,
                                                                         Vector3D<Real_v> const &point,
                                                                         Vector3D<Real_v> const &direction, Real_v step,
                                                                         LeafFunction &&visit)
  {
    unsigned int stack[FacetBVH_t::kTraversalStackSize];
    int stackSize = 0;
    const Vector3D<float> bvhPoint(static_cast<float>(point[0]), static_cast<float>(point[1]),
                                   static_cast<float>(point[2]));
    const Vector3D<float> invDirection(1.f / vecgeom::NonZero(direction[0]), 1.f / vecgeom::NonZero(direction[1]),
                                       1.f / vecgeom::NonZero(direction[2]));
    const unsigned int octant =
        (invDirection[0] < 0.f ? 1u : 0u) | (invDirection[1] < 0.f ? 2u : 0u) | (invDirection[2] < 0.f ? 4u : 0u);
    float bvhStep = static_cast<float>(step);
    float approach;
    if (!bvh.IntersectNode(0, bvhPoint, invDirection, octant, bvhStep, approach)) return;

    auto const *nodes  = bvh.GetNodes();
    auto const *primId = bvh.GetPrimId();
    unsigned int top   = 0;
    for (;;) {
      while (!nodes[top].IsLeaf()) {
        const unsigned int childL = nodes[top].first;
        const unsigned int childR = childL + 1;
        float approachL, approachR;
        const bool hitL = bvh.IntersectNode(childL, bvhPoint, invDirection, octant, bvhStep, approachL);
        const bool hitR = bvh.IntersectNode(childR, bvhPoint, invDirection, octant, bvhStep, approachR);
        if (hitL) {
          unsigned int near = childL;
          if (hitR) {
            unsigned int far = childR;
            if (approachR < approachL) {
              near = childR;
              far  = childL;
            }
            stack[stackSize++] = far;
          }
          top = near;
        } else if (hitR) {
          top = childR;
        } else {
          goto pop;
        }
      }

      for (int i = 0; i < nodes[top].count; ++i) {
        const int facetId = primId[nodes[top].first + i];
        approach          = 0.f;
        if (check_leaf_bb && !bvh.IntersectBox(facetId, bvhPoint, invDirection, bvhStep, approach)) continue;
        if (visit(facetId, bvhStep, approach)) {
          if constexpr (stop_at_hit) return;
          break;
        }
      }

    pop:
      if (stackSize == 0) return;
      top = stack[--stackSize];
    }
  }

  /**
   * @brief Traverse facet boxes with a mutable squared-Safety bound.
   * @details Children are visited near-first. Entry queries may accept a
   * conservative node estimate instead of descending to facets.
   * @tparam return_top_estimate Permit conservative node estimates.
   * @tparam LeafFunction Facet Safety visitor type.
   * @param bvh Facet hierarchy.
   * @param localpoint Local query point.
   * @param limitSq Initial squared-Safety bound.
   * @param visit Exact facet Safety operation receiving `(facetId, safetySqr)`;
   *              it may reduce the squared bound by reference.
   * @return Reduced or conservatively estimated squared Safety.
   */
  template <bool return_top_estimate, typename LeafFunction>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static Precision VisitPointSafetyBoxes(
      FacetBVH_t const &bvh, Vector3D<Precision> const &localpoint, Precision limitSq, LeafFunction &&visit)
  {
    struct StackItem {
      unsigned int node;
      float safetySqr;
    };
    StackItem stack[FacetBVH_t::kTraversalStackSize], *stackPtr = &stack[1];
    const Vector3D<float> point(static_cast<float>(localpoint[0]), static_cast<float>(localpoint[1]),
                                static_cast<float>(localpoint[2]));
    float safetySqr = static_cast<float>(limitSq);
    stack[0]        = StackItem{0, bvh.NodeSafetySqr(0, point)};
    if constexpr (return_top_estimate) {
      if (stack[0].safetySqr > 0.f) return static_cast<Precision>(stack[0].safetySqr);
    }

    auto const *nodes  = bvh.GetNodes();
    auto const *primId = bvh.GetPrimId();
    do {
      const auto current     = *--stackPtr;
      const unsigned int top = current.node;
      if (current.safetySqr > safetySqr) continue;
      if (nodes[top].IsLeaf()) {
        for (int i = 0; i < nodes[top].count; ++i)
          visit(primId[nodes[top].first + i], safetySqr);
        continue;
      }

      const unsigned int childL = nodes[top].first;
      const unsigned int childR = childL + 1;
      const float safetyL       = bvh.NodeSafetySqr(childL, point);
      const float safetyR       = bvh.NodeSafetySqr(childR, point);
      bool traverseL            = safetyL < safetySqr;
      bool traverseR            = safetyR < safetySqr;
      if constexpr (return_top_estimate) {
        constexpr float threshold = 0.1f;
        if (traverseL && safetyL > threshold) {
          traverseL = false;
          safetySqr = vecCore::math::Min(safetySqr, safetyL);
        }
        if (traverseR && safetyR > threshold) {
          traverseR = false;
          safetySqr = vecCore::math::Min(safetySqr, safetyR);
        }
      }
      if (safetyR < safetyL) {
        if (traverseL) *stackPtr++ = StackItem{childL, safetyL};
        if (traverseR) *stackPtr++ = StackItem{childR, safetyR};
      } else {
        if (traverseR) *stackPtr++ = StackItem{childR, safetyR};
        if (traverseL) *stackPtr++ = StackItem{childL, safetyL};
      }
    } while (stackPtr > stack);
    return FacetBVH_t::FinalizeSafetySqr(safetySqr);
  }

  // Repeat only an ambiguous parity query with consistent half-open ownership.
  // Inline the Boolean fallback so register lifetimes can be reused between
  // the ordinary traversal and this rarely executed ownership calculation.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static bool ResolveParity(UnplacedStruct_t const &tessellated,
                                                                         Vector3D<Real_v> const &point, Real_v rayEPS)
  {
    int count = 0;
    VisitRayBoxes<kCheckFacetBB>(
        *tessellated.fBVH, point, tessellated.fTestDir, InfinityLength<Real_v>(), [&](int primID, float &, float) {
          if (tessellated.fFacets[primID].OwnsParityCrossing(point, tessellated.fTestDir, rayEPS)) ++count;
          return false;
        });
    return count % 2 == 1;
  }

  /// @brief Test whether a point is inside the closed tessellated shell.
  /// @details Uses the shared first-crossing classifier without a surface-band
  /// exclusion. Retains the original extent and positive crossing cutoff.
  /// @tparam Real_v Floating-point scalar type.
  /// @param tessellated Runtime tessellated data.
  /// @param point Local point to test.
  /// @param[out] contains Containment decision, including accepted surface-band points.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void Contains(UnplacedStruct_t const &tessellated,
                                                                    Vector3D<Real_v> const &point, bool &contains)
  {
    bool unused = false, outside = false;
    GenericKernelForContainsAndInside<Real_v, false>(tessellated, point, unused, outside);
    contains = !outside;
  }

  /// @brief Test the spatial surface band of a finite triangle.
  /// @details Keep this check out of line on CUDA for the measured distance-query benefit.
  /// @tparam check_edge_band Also resolve exterior projections near finite edges.
  /// @param facet Triangle with a unit normal.
  /// @param point Local point to classify.
  /// @return True for a certified point in the finite triangle's half-tolerance band.
  template <typename Real_v, bool check_edge_band = true>
#ifdef VECCORE_CUDA_DEVICE_COMPILATION
  VECCORE_FORCE_NOINLINE
#else
  VECGEOM_FORCE_INLINE
#endif
      VECCORE_ATT_HOST_DEVICE static bool IsOnFiniteFacet(TriangularTile<Precision> const &facet,
                                                          Vector3D<Real_v> const &point)
  {
    const auto planeDistance = (point - facet.fVertices[0]).Dot(facet.fNormal);
    if (Abs(planeDistance) >= kHalfTolerance) return false;
    for (int i = 0; i < 3; ++i) {
      const auto edge = facet.fVertices[(i + 1) % 3] - facet.fVertices[i];
      const auto side = edge.Cross(facet.fNormal).Dot(point - facet.fVertices[i]);
      if (side > Real_v(0.)) {
        if constexpr (!check_edge_band) return false;
        if (side * side >= Real_v(kHalfTolerance * kHalfTolerance) * edge.Mag2()) return false;
        // An exterior projection is nearest to the boundary. Expanded edge
        // planes overestimate its band at acute vertices; check finite edges.
        for (int j = 0; j < 3; ++j) {
          if (DistanceToLineSegmentSquared<Real_v>(facet.fVertices[j], facet.fVertices[(j + 1) % 3], point) <
              Real_v(kHalfTolerance * kHalfTolerance))
            return true;
        }
        return false;
      }
    }
    return true;
  }

  /// @brief Classify a point by its first oriented ray crossing.
  /// @details Any finite-facet surface hit completes classification immediately.
  /// Otherwise the first unambiguous outward crossing starts inside; an inward
  /// crossing or no crossing starts outside. Uncertain edge ownership, tangency
  /// or tied opposite crossings fall back to the full parity classifier.
  /// @tparam Real_v Floating-point scalar type.
  /// @param tessellated Runtime tessellated data.
  /// @param point Local point to classify.
  /// @param[out] inside Set to `kInside`, `kOutside`, or `kSurface`.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void Inside(UnplacedStruct_t const &tessellated,
                                                                  Vector3D<Real_v> const &point, Inside_t &inside)
  {
    bool completelyinside = false, completelyoutside = false;
    GenericKernelForContainsAndInside<Real_v, true>(tessellated, point, completelyinside, completelyoutside);
    inside = completelyoutside ? kOutside : (completelyinside ? kInside : kSurface);
  }

  /// @brief Shared extent rejection and first-crossing classification.
  /// @tparam Real_v Floating-point scalar type.
  /// @tparam ForInside Resolve the surface band; otherwise retain Contains' crossing convention.
  /// @param tessellated Runtime tessellated data.
  /// @param point Local query point.
  /// @param[out] completelyinside True for a classified interior point.
  /// @param[out] completelyoutside True when the selected query rejects the point.
  template <typename Real_v, bool ForInside>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void GenericKernelForContainsAndInside(
      UnplacedStruct_t const &tessellated, Vector3D<Real_v> const &point, bool &completelyinside,
      bool &completelyoutside)
  {
    constexpr Real_v margin = ForInside ? Real_v(kHalfTolerance) : Real_v(0.);
    bool contains           = false;
    ABBoxImplementation::ABBoxContainsKernel(tessellated.fMinExtent - Vector3D<Real_v>(margin),
                                             tessellated.fMaxExtent + Vector3D<Real_v>(margin), point, contains);
    completelyinside  = false;
    completelyoutside = !contains;
    if (!contains) return;

    Inside_t inside;
    ClassifyInside<Real_v, ForInside>(tessellated, point, inside);
    completelyinside  = inside == kInside;
    completelyoutside = inside == kOutside;
  }

  /// @brief Classify a point accepted by the extent check.
  /// @details Separate traversal state from the bounding-box fast path and let
  /// the compiler decide whether to inline the shared classifier.
  /// @param tessellated Runtime tessellated data.
  /// @tparam ForInside Enable the finite surface band and Inside crossing cutoff.
  /// @param point Local point accepted by the selected extent check.
  /// @param[out] inside Classification; Contains mode produces only kInside or kOutside.
  template <typename Real_v, bool ForInside>
  VECCORE_ATT_HOST_DEVICE static void ClassifyInside(UnplacedStruct_t const &tessellated, Vector3D<Real_v> const &point,
                                                     Inside_t &inside)
  {
    Real_v nearest     = InfinityLength<Real_v>();
    Real_v orientation = 0.;
    bool onSurface     = false;
    bool ambiguous     = false;
    const auto &bounds = tessellated.fBVH->GetNodes()[0].bounds;
    const float extent =
        Max(bounds[1] - bounds[0], Max(bounds[3] - bounds[2], bounds[5] - bounds[4])) + float(kTolerance);
    // Same triple-product error budget as Contains, with the surface allowance.
    const float parityError =
        256.f * float(kEpsilonT<Real_v>) * extent * extent * Vector3D<float>(tessellated.fTestDir).Abs().Max();
    // Contains preserves its original positive ray cutoff and never removes
    // the inner surface band merely because Inside would return kSurface.
    constexpr Real_v rayEPS = ForInside ? Real_v(-kHalfTolerance) : Real_v(1.e-8);
    float surfaceApproach   = 0.f;
    if constexpr (ForInside) {
      const Vector3D<float> bvhPoint(point);
      const Vector3D<float> invDirection(1.f / NonZero(tessellated.fTestDir[0]), 1.f / NonZero(tessellated.fTestDir[1]),
                                         1.f / NonZero(tessellated.fTestDir[2]));
      // A surface point must also be close to the enclosing facet box. Include
      // the point's float-rounding error in this rejection-only bound; retain
      // the exact finite-triangle test and its original tolerance below.
      const float surfaceProbe = Max(float(kTolerance), 4.f * kEpsilonT<float> * bvhPoint.Abs().Max());
      // A point within surfaceProbe of a box cannot approach it later than
      // surfaceProbe / min(abs(direction)). Near-zero components only loosen
      // this rejection, so grazing test rays retain the exact surface check.
      surfaceApproach = surfaceProbe * invDirection.Abs().Max();
    }
    VisitRayBoxes<false, true>(*tessellated.fBVH, point, tessellated.fTestDir, InfinityLength<Real_v>(),
                               [&](int primID, float &bvhStep, float approach) {
                                 const auto &facet = tessellated.fFacets[primID];
                                 if constexpr (ForInside) {
                                   if (approach < surfaceApproach && IsOnFiniteFacet(facet, point)) {
                                     onSurface = true;
                                     return true;
                                   }
                                 }
                                 Real_v distance;
                                 if (facet.template ParityCrossing<true>(point, tessellated.fTestDir, rayEPS, ambiguous,
                                                                         parityError, &distance)) {
                                   const Real_v crossing = facet.fNormal.Dot(tessellated.fTestDir);
                                   if (Abs(distance - nearest) < kTolerance && crossing * orientation < 0.)
                                     ambiguous = true;
                                   if (distance < nearest) {
                                     nearest     = distance;
                                     orientation = crossing;
                                     // Retain nearby boxes needed for the independent surface test,
                                     // as well as tied crossings and float traversal roundoff.
                                     bvhStep = Max(surfaceApproach,
                                                   static_cast<float>(Max(distance, Real_v(0.)) + kTolerance) *
                                                       (1.f + 16.f * kEpsilonT<float>));
                                   }
                                 }
                                 return false;
                               });
    if constexpr (ForInside) {
      if (onSurface) {
        inside = kSurface;
        return;
      }
    }
    if (ambiguous ||
        (nearest < InfinityLength<Real_v>() && (nearest < kTolerance || Abs(orientation) < 32. * kEpsilonT<Real_v>))) {
      // The traversal retained every surface-probe box. Only crossing ownership
      // remains unresolved, so reuse the existing half-open parity pass.
      inside = ResolveParity(tessellated, point, rayEPS) ? Inside_t(kInside) : Inside_t(kOutside);
      return;
    }
    inside = orientation > 0. ? Inside_t(kInside) : Inside_t(kOutside);
  }

  /// @brief Bound roundoff in a ray's barycentric triple products.
  /// @details L bounds vertex-vertex and vertex-point component differences.
  /// For unit navigation directions, 6*gamma(8)*L^2 bounds each triple
  /// product; 256*epsilon also covers differences of three products and the
  /// float evaluation of this rejection-only bound. Inside/surface starts use
  /// the outward-rounded root; entry queries additionally include their origin.
  /// @param tessellated Mesh and its outward-rounded BVH bounds.
  /// @param point Ray origin.
  /// @param fromInside Whether the origin is already known to be inside/on the mesh.
  /// @return Conservative barycentric uncertainty, not a spatial tolerance.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static Real_v RayIntersectionError(UnplacedStruct_t const &tessellated,
                                                                                  Vector3D<Real_v> const &point,
                                                                                  bool fromInside = false)
  {
    const auto &b = tessellated.fBVH->GetNodes()[0].bounds;
    float length  = Max(b[1] - b[0], Max(b[3] - b[2], b[5] - b[4])) + float(kTolerance);
    if (!fromInside) {
      const Vector3D<float> p(point), lo(b[0], b[2], b[4]), hi(b[1], b[3], b[5]);
      length = Max(length, Max((p - lo).Max(), (hi - p).Max()));
    }
    return 256.f * float(kEpsilonT<Real_v>) * length * length;
  }

  /// @brief Compute entry distance without wrong-side classification.
  /// @details The BVH traversal rejects facets whose normal faces away from an
  /// incoming ray. Surface starts can be recorded as a zero-entry fallback, but
  /// an ordinary positive intersection within tolerance wins when present.
  /// @tparam Real_v Floating-point scalar type.
  /// @param tessellated Runtime tessellated data.
  /// @param point Local start point.
  /// @param direction Unit local direction.
  /// @param stepMax Maximum distance to consider.
  /// @param[out] distance Nearest entry candidate or `kInfLength`.
  /// @param allow_surface_zero Allow a tolerated surface start to return zero.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void DistanceToInNoConv(UnplacedStruct_t const &tessellated,
                                                                              Vector3D<Real_v> const &point,
                                                                              Vector3D<Real_v> const &direction,
                                                                              Real_v const &stepMax, Real_v &distance,
                                                                              bool allow_surface_zero = false)
  {
    distance                        = InfinityLength<Real_v>();
    bool found_surface_entry        = false;
    const auto intersectionError    = RayIntersectionError(tessellated, point);
    Real_v surface_entry_projection = Real_v(0.);

    // NOTE: a quick intersection check against the outer bounding box is already done as part of the BVH
    // intersection and does not need to be done in addition here
    VisitRayBoxes<kCheckFacetBB>(*tessellated.fBVH, point, direction, stepMax, [&](int primID, float &bvhStep, float) {
      const auto &facet = tessellated.fFacets[primID];
      // we are checking a triangle. Rule out early by a simple normal check
      const auto sp = facet.fNormal.Dot(direction);
      if (allow_surface_zero && sp < -kToleranceDist<Real_v>) {
        const auto plane_dist = (point - facet.fVertices[0]).Dot(facet.fNormal);
        if (Abs(plane_dist) < kToleranceDist<Real_v> &&
            facet.template SafetySq<Real_v>(point) < kToleranceDistSquared<Real_v>) {
          // Keep this as a fallback only: a positive ordinary intersection,
          // when found within normal tolerance, must win over the surface-zero
          // convention.
          found_surface_entry      = true;
          surface_entry_projection = Max(surface_entry_projection, -sp);
        }
      }
      const bool wrong_orientation = sp > 0.; // coming from outside the dot product must be negative
      if (wrong_orientation) return false;
      // Match the exit tangent convention at a surface start. The triangle
      // determinant scales with facet area and can admit an unresolved normal
      // crossing. Exclude its zero entry with the existing root cutoff, while
      // keeping later crossings eligible and avoiding a second hit filter.
      const Real_v entryMin =
          allow_surface_zero && !(sp < -kToleranceDist<Real_v>) ? kToleranceDist<Real_v> : -Real_v(kHalfTolerance);
      const auto this_dist = facet.template Distance<true, 1>(point, direction, entryMin, intersectionError);
      if (this_dist < distance) {
        // update stuff
        distance = Max(this_dist, Real_v(0.));
        bvhStep  = static_cast<float>(this_dist); // prune later boxes in BVH precision
      }
      return false;
    });
    if (found_surface_entry &&
        (distance == InfinityLength<Real_v>() || distance * surface_entry_projection > kToleranceDist<Real_v>)) {
      distance = Real_v(0.);
    }
  }

  /// @brief Compute distance from an outside or surface point to enter.
  /// @details Performs public wrong-side classification before delegating to
  /// the no-convention helper. Inside starts return `-1`.
  /// @note Known limitation: a ray can report a finite edge/vertex contact
  /// while remaining outside on both sides. Entry does not apply the local
  /// material-continuation test used to resolve concave exit touches.
  /// @tparam Real_v Floating-point scalar type.
  /// @param tessellated Runtime tessellated data.
  /// @param point Local start point.
  /// @param direction Unit local direction.
  /// @param stepMax Maximum distance to consider.
  /// @param[out] distance Entry distance, `-1`, or `kInfLength`.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void DistanceToIn(UnplacedStruct_t const &tessellated,
                                                                        Vector3D<Real_v> const &point,
                                                                        Vector3D<Real_v> const &direction,
                                                                        Real_v const &stepMax, Real_v &distance)
  {
#ifndef VECGEOM_NO_WRONGSIDE_CONV
    Inside_t inside;
    Inside<Real_v>(tessellated, point, inside);
    if (inside == kInside) {
      distance = Real_v(-1.);
      return;
    }
    const bool allow_surface_zero = inside == kSurface;
#else
    const bool allow_surface_zero = false;
#endif
    DistanceToInNoConv<Real_v>(tessellated, point, direction, stepMax, distance, allow_surface_zero);
  }

  /// @brief Test whether a triangle hit is within tolerance of an edge line.
  /// @param facet Triangle already selected by the distance reduction.
  /// @param hit Point on that triangle, within distance-query roundoff.
  /// @return Whether local edge/vertex qualification is required.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static bool IsNearFacetEdge(TriangularTile<Precision> const &facet,
                                                                           Vector3D<Real_v> const &hit)
  {
    bool near = false;
    for (int i = 0; i < 3; ++i) {
      const auto edge = facet.fVertices[(i + 1) % 3] - facet.fVertices[i];
      const auto side = edge.Cross(facet.fNormal).Dot(hit - facet.fVertices[i]);
      near |= side * side < kToleranceDistSquared<Real_v> * edge.Mag2();
    }
    return near;
  }

  /// @brief Test incidence using plane distance and finite-triangle safety.
  /// @param facet Candidate triangle.
  /// @param hit Candidate boundary point.
  /// @return Whether the point is within the triangle's spatial tolerance.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static bool IsIncidentFacet(TriangularTile<Precision> const &facet,
                                                                           Vector3D<Real_v> const &hit)
  {
    return Abs((hit - facet.fVertices[0]).Dot(facet.fNormal)) < kToleranceDist<Real_v> &&
           facet.template SafetySq<Real_v>(hit) < kToleranceDistSquared<Real_v>;
  }

  /// @brief Compute a corner's angle for the local angle-weighted normal.
  /// @param facet Incident triangle.
  /// @param vertex Corner index.
  /// @return Interior angle in radians.
  /// @details Isolate the rare angle calculation to limit live register state.
  template <typename Real_v>
  VECCORE_FORCE_NOINLINE VECCORE_ATT_HOST_DEVICE static Real_v FacetVertexAngle(TriangularTile<Precision> const &facet,
                                                                                int vertex)
  {
    const auto a = facet.fVertices[(vertex + 1) % 3] - facet.fVertices[vertex];
    const auto b = facet.fVertices[(vertex + 2) % 3] - facet.fVertices[vertex];
    return ATan2(a.Cross(b).Mag(), a.Dot(b));
  }

  /// @brief Find the direction's nearest feature on an incident tangent cone.
  /// @param facet Incident triangle.
  /// @param hit Boundary point determining a full plane, half-plane, or vertex wedge.
  /// @param direction Outgoing unit direction from the boundary point.
  /// @param[out] distanceSq Squared distance from direction to that cone.
  /// @param[out] signedProjection Angle-weighted outward-normal projection.
  /// @details Project onto the plane, edge line, or vertex half-rays. Weight
  /// the nearest feature by 2*pi (face), pi (edge), or its corner angle (vertex).
  /// No displaced position or finite-distance classification is evaluated.
  /// Separating face/edge and vertex arithmetic bounds the live register state.
  template <typename Real_v>
  VECCORE_FORCE_NOINLINE VECCORE_ATT_HOST_DEVICE static void FacetConeDistance(TriangularTile<Precision> const &facet,
                                                                               Vector3D<Real_v> const &hit,
                                                                               Vector3D<Real_v> const &direction,
                                                                               Real_v &distanceSq,
                                                                               Real_v &signedProjection)
  {
    const Real_v roundoff = 64. * kEpsilonT<Real_v>;
    int vertex            = -1;
    const auto normalPart = facet.fNormal.Dot(direction);
    for (int i = 0; i < 3; ++i) {
      if ((hit - facet.fVertices[i]).Mag2() < kToleranceDistSquared<Real_v>) {
        vertex = i;
        break;
      }
    }
    distanceSq    = normalPart * normalPart;
    Real_v weight = Real_v(kTwoPi);
    if (vertex < 0) {
      for (int i = 0; i < 3; ++i) {
        const auto edge       = facet.fVertices[(i + 1) % 3] - facet.fVertices[i];
        const auto sideVector = edge.Cross(facet.fNormal);
        const auto side       = sideVector.Dot(hit - facet.fVertices[i]);
        if (side * side < kToleranceDistSquared<Real_v> * edge.Mag2()) {
          const auto motion = sideVector.Dot(direction);
          const auto error  = roundoff * edge.Abs().Max();
          if (motion > error) distanceSq = (direction - (direction.Dot(edge) / edge.Mag2()) * edge).Mag2();
          if (motion >= -error) weight = Real_v(kPi);
          break;
        }
      }
    } else {
      const auto a       = facet.fVertices[(vertex + 1) % 3] - facet.fVertices[vertex];
      const auto b       = facet.fVertices[(vertex + 2) % 3] - facet.fVertices[vertex];
      const auto motionA = a.Cross(facet.fNormal).Dot(direction);
      const auto motionB = -b.Cross(facet.fNormal).Dot(direction);
      const auto errorA = roundoff * a.Abs().Max(), errorB = roundoff * b.Abs().Max();
      bool atVertex;
      if (motionA > errorA || motionB > errorB) {
        const auto pa = Max(Real_v(0.), direction.Dot(a) / a.Mag2());
        const auto pb = Max(Real_v(0.), direction.Dot(b) / b.Mag2());
        distanceSq    = Min((direction - pa * a).Mag2(), (direction - pb * b).Mag2());
        atVertex      = pa == Real_v(0.) && pb == Real_v(0.);
        weight        = Real_v(kPi);
      } else {
        atVertex = (direction - normalPart * facet.fNormal).Mag2() <= roundoff * roundoff * direction.Mag2();
        if (Abs(motionA) <= errorA || Abs(motionB) <= errorB) weight = Real_v(kPi);
      }
      if (atVertex) weight = FacetVertexAngle<Real_v>(facet, vertex);
    }
    signedProjection = weight * normalPart;
  }

  /// @brief Test whether a boundary direction continues locally in material.
  /// @param tessellated Closed, consistently oriented mesh.
  /// @param hit Non-supporting edge/vertex hit.
  /// @param direction Outgoing unit direction.
  /// @return True when the nearest tangent-cone features have inward signed distance.
  /// @details Query only incident facets and combine their angle-weighted
  /// normals at tied nearest features. This is a local tangent-cone predicate,
  /// not an Inside query or displaced-point probe; tangent directions return false.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static bool ForwardInside(UnplacedStruct_t const &tessellated,
                                                                         Vector3D<Real_v> const &hit,
                                                                         Vector3D<Real_v> const &direction)
  {
    Real_v closest = InfinityLength<Real_v>(), signedProjection = 0.;
    const Real_v roundoff = 64. * kEpsilonT<Real_v>;
    VisitPointSafetyBoxes<false>(*tessellated.fBVH, hit, kToleranceDistSquared<Real_v>, [&](int primID, float &) {
      const auto &facet = tessellated.fFacets[primID];
      if (!IsIncidentFacet(facet, hit)) return;
      Real_v distanceSq, projection;
      FacetConeDistance(facet, hit, direction, distanceSq, projection);
      const auto tie = roundoff * Max(distanceSq, closest);
      if (distanceSq < closest - tie) {
        closest          = distanceSq;
        signedProjection = projection;
      } else if (Abs(distanceSq - closest) <= tie) {
        signedProjection += projection;
      }
    });
    return signedProjection < -roundoff;
  }

  /// @brief Find the next outward facet beyond the current incident feature.
  /// @param tessellated Mesh and BVH.
  /// @param hit Current boundary point.
  /// @param direction Unit ray direction.
  /// @param[out] next Distance beyond hit, or infinity when none is found.
  /// @param[out] nextFacet Selected facet index, or -1 for a miss.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void NextNonincidentExit(UnplacedStruct_t const &tessellated,
                                                                               Vector3D<Real_v> const &hit,
                                                                               Vector3D<Real_v> const &direction,
                                                                               Real_v &next, int &nextFacet)
  {
    next                         = InfinityLength<Real_v>();
    nextFacet                    = -1;
    const auto intersectionError = RayIntersectionError(tessellated, hit, true);
    VisitRayBoxes<kCheckFacetBB>(*tessellated.fBVH, hit, direction, next, [&](int primID, float &bvhStep, float) {
      const auto &other = tessellated.fFacets[primID];
      if (other.fNormal.Dot(direction) < 0.) return false;
      // Skip the current incident root, including shallow-ray roundoff.
      if (IsIncidentFacet(other, hit)) return false;
      const auto candidate = other.template Distance<true, -1>(hit, direction, Real_v(0.), intersectionError);
      if (candidate < next) {
        next      = candidate;
        nextFacet = primID;
        bvhStep   = static_cast<float>(candidate);
      }
      return false;
    });
  }

  /// @brief Skip concave touches until the ray leaves material.
  /// @param tessellated Mesh and BVH.
  /// @param hit Current candidate hit.
  /// @param direction Unit ray direction.
  /// @param[in,out] distance Accumulated distance to the selected exit.
  /// @param[in,out] facetIndex Selected exit facet, or -1 if no exit is found.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void ResolveEdgeExit(UnplacedStruct_t const &tessellated,
                                                                           Vector3D<Real_v> hit,
                                                                           Vector3D<Real_v> const &direction,
                                                                           Real_v &distance, int &facetIndex)
  {
    while (ForwardInside(tessellated, hit, direction)) {
      Real_v next;
      int nextFacet;
      NextNonincidentExit(tessellated, hit, direction, next, nextFacet);
      if (nextFacet < 0) {
        distance   = InfinityLength<Real_v>();
        facetIndex = -1;
        return;
      }
      distance += next;
      facetIndex = nextFacet;
      hit += next * direction;
      if (tessellated.FacetSupports(facetIndex) || !IsNearFacetEdge(tessellated.fFacets[facetIndex], hit)) return;
    }
  }

  /// @brief Compute exit distance without wrong-side classification.
  /// @details The BVH traversal considers outward-facing facets and accepts
  /// tolerated tangent/outward starts as zero exits on the current facet.
  /// @tparam Real_v Floating-point scalar type.
  /// @param tessellated Runtime tessellated data.
  /// @param point Local start point.
  /// @param direction Unit local direction.
  /// @param stepMax Maximum distance to consider.
  /// @param[out] distance Nearest exit candidate, `stepMax` when no exit is
  /// found before a finite limit, or `kInfLength` for an unbounded miss.
  /// @param[out] hitFacet Selected facet index, or -1 if no candidate is accepted;
  /// required only when RecordFacet is true.
  template <typename Real_v, bool AllowSurfaceZero, bool RecordFacet = false>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void DistanceToOutNoConv(UnplacedStruct_t const &tessellated,
                                                                               Vector3D<Real_v> const &point,
                                                                               Vector3D<Real_v> const &direction,
                                                                               Real_v const &stepMax, Real_v &distance,
                                                                               int *hitFacet = nullptr)
  {
    distance                     = InfinityLength<Real_v>();
    int nearestFacet             = -1;
    const auto intersectionError = RayIntersectionError(tessellated, point, true);
    if constexpr (RecordFacet) *hitFacet = -1;

    VisitRayBoxes<kCheckFacetBB>(*tessellated.fBVH, point, direction, stepMax, [&](int primID, float &bvhStep, float) {
      const auto &facet = tessellated.fFacets[primID];
      // we are checking a triangle. Rule out early by a simple normal check
      const auto sp                = facet.fNormal.Dot(direction);
      const bool wrong_orientation = sp < 0.; // coming from inside the dot product must be positive
      const auto plane_dist        = (point - facet.fVertices[0]).Dot(facet.fNormal);
      // Accept t=0 for tangent/outward surface starts on the current facet;
      // genuinely inward starts must continue to the next exit surface.
      if (AllowSurfaceZero && sp >= -kToleranceDist<Real_v> && Abs(plane_dist) < kToleranceDist<Real_v>) {
        const auto safety_sq = facet.template SafetySq<Real_v>(point);
        if (safety_sq < kToleranceDistSquared<Real_v>) {
          distance     = Real_v(0.);
          bvhStep      = 0.f;
          nearestFacet = primID;
          return true;
        }
      }
      if (wrong_orientation) return false;
      // the -kHalfTolerance is to get 0 if we are on surface
      const auto this_dist = facet.template Distance<true, -1>(point, direction, -kHalfTolerance, intersectionError);
      if (this_dist < distance) {
        // update stuff
        distance     = Max(this_dist, Real_v(0.));
        bvhStep      = static_cast<float>(this_dist); // prune later boxes in BVH precision
        nearestFacet = primID;
      }
      return false;
    });
    if (nearestFacet >= 0 && !tessellated.FacetSupports(nearestFacet)) {
      const auto hit = point + distance * direction;
      if (IsNearFacetEdge(tessellated.fFacets[nearestFacet], hit)) {
        const auto originalDistance = distance;
        ResolveEdgeExit(tessellated, hit, direction, distance, nearestFacet);
        if (distance > stepMax && distance > originalDistance) {
          distance     = stepMax;
          nearestFacet = -1;
        }
      }
    }
    if constexpr (RecordFacet) *hitFacet = nearestFacet;
    if (distance == InfinityLength<Real_v>() && stepMax < InfinityLength<Real_v>()) {
      distance = stepMax;
    }
  }

  /// @brief Compute distance from an inside or surface point to leave.
  /// @details Performs public wrong-side classification before delegating to
  /// the no-convention helper. Outside starts return `-1`.
  /// @tparam Real_v Floating-point scalar type.
  /// @param tessellated Runtime tessellated data.
  /// @param point Local start point.
  /// @param direction Unit local direction.
  /// @param stepMax Maximum distance to consider.
  /// @param[out] distance Exit distance, `-1`, `stepMax` for finite-limit
  /// misses, or `kInfLength`.
  /// @param[out] hitInfo Optional supporting exit normal; false validity leaves normal storage untouched.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void DistanceToOut(UnplacedStruct_t const &tessellated,
                                                                         Vector3D<Real_v> const &point,
                                                                         Vector3D<Real_v> const &direction,
                                                                         Real_v const &stepMax, Real_v &distance,
                                                                         SurfaceHitView<Real_v> *hitInfo = nullptr)
  {
    if (hitInfo) {
      hitInfo->Clear();
      int hitFacet;
      DistanceToOutWithFacet(tessellated, point, direction, stepMax, distance, hitFacet);
      if (hitFacet >= 0) {
        const bool leftBehind = tessellated.FacetSupports(hitFacet);
        hitInfo->SetLeftBehind(leftBehind);
        if (leftBehind && hitInfo->WantsNormal())
          hitInfo->SetNormal(Vector3D<Real_v>(tessellated.fFacets[hitFacet].fNormal));
      }
      return;
    }
#ifndef VECGEOM_NO_WRONGSIDE_CONV
    Inside_t inside;
    Inside<Real_v>(tessellated, point, inside);
    if (inside == kOutside) {
      distance = Real_v(-1.);
      return;
    }
    if (inside == kSurface) {
      DistanceToOutNoConv<Real_v, true>(tessellated, point, direction, stepMax, distance);
      return;
    }
#endif
    DistanceToOutNoConv<Real_v, false>(tessellated, point, direction, stepMax, distance);
  }

  /// @brief Compute an exit distance and return the facet selected by the existing BVH minimization.
  /// @param tessellated Compact mesh data and facet BVH.
  /// @param point Local start point.
  /// @param direction Unit local direction.
  /// @param stepMax Maximum distance to consider.
  /// @param[out] distance Exit distance with the ordinary wrong-side convention.
  /// @param[out] hitFacet Winning facet, or -1 for an unaccepted crossing.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void DistanceToOutWithFacet(UnplacedStruct_t const &tessellated,
                                                                                  Vector3D<Real_v> const &point,
                                                                                  Vector3D<Real_v> const &direction,
                                                                                  Real_v const &stepMax,
                                                                                  Real_v &distance, int &hitFacet)
  {
    hitFacet = -1;
#ifndef VECGEOM_NO_WRONGSIDE_CONV
    Inside_t inside;
    Inside<Real_v>(tessellated, point, inside);
    if (inside == kOutside) {
      distance = Real_v(-1.);
      return;
    }
    if (inside == kSurface) {
      DistanceToOutNoConv<Real_v, true, true>(tessellated, point, direction, stepMax, distance, &hitFacet);
      return;
    }
#endif
    DistanceToOutNoConv<Real_v, false, true>(tessellated, point, direction, stepMax, distance, &hitFacet);
  }

  /// @brief Compute unsigned safety to enter without wrong-side classification.
  /// @details Uses cached surface anchor points to seed an upper BVH query
  /// limit, then asks the facet BVH for a squared safety estimate.
  /// @tparam Real_v Floating-point scalar type.
  /// @param tessellated Runtime tessellated data.
  /// @param point Local point.
  /// @param[out] safety Non-negative entry safety.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SafetyToInNoConv(UnplacedStruct_t const &tessellated,
                                                                            Vector3D<Real_v> const &point,
                                                                            Real_v &safety)
  {
    int isurf;

    // get a quick upper limit from the min-distance to fixed set of anchor points on the surface
    // --> this limits BVH search from the start
    // TODO: this should also be vectorizable on the CPU
    float upper_limit_sq = InfinityLength<float>();
    const float px = point.x(), py = point.y(), pz = point.z();
    for (int i = 0; i < TessellatedRuntimeStruct<float>::N; ++i) {
      float dx       = tessellated.fTestPoints_x[i] - px;
      float dy       = tessellated.fTestPoints_y[i] - py;
      float dz       = tessellated.fTestPoints_z[i] - pz;
      float dist2    = dx * dx + dy * dy + dz * dz;
      upper_limit_sq = Min(dist2, upper_limit_sq);
    }

    constexpr bool approxSafety = true; // can return early (without detailed safety, but never 0)
    const Real_v safetysq       = SafetySq<Real_v, approxSafety, double>(tessellated, point, isurf, upper_limit_sq);
    safety                      = Sqrt(safetysq);

    // Keep the NoConv helper unsigned; public wrappers own wrong-side conventions.
    // safety on boundary should be zero --> see ShapeTester
    if (safety < kTolerance) {
      safety = 0.;
    }
  }

  /// @brief Compute safety from an outside or surface point to enter.
  /// @details Performs public wrong-side classification before delegating to
  /// the no-convention helper. Inside starts return `-1`, and surface starts
  /// return zero.
  /// @tparam Real_v Floating-point scalar type.
  /// @param tessellated Runtime tessellated data.
  /// @param point Local point.
  /// @param[out] safety Entry safety, `-1`, or zero for surface starts.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SafetyToIn(UnplacedStruct_t const &tessellated,
                                                                      Vector3D<Real_v> const &point, Real_v &safety)
  {
#ifndef VECGEOM_NO_WRONGSIDE_CONV
    Inside_t inside;
    Inside<Real_v>(tessellated, point, inside);
    if (inside == kInside) {
      safety = Real_v(-1.);
      return;
    }
    // Surface-point safety is defined as zero; avoid the expensive BVH safety query in this case.
    if (inside == kSurface) {
      safety = Real_v(0.);
      return;
    }
#endif
    SafetyToInNoConv<Real_v>(tessellated, point, safety);
  }

  /// @brief Compute unsigned safety to leave without wrong-side classification.
  /// @details Uses cached surface anchor points to seed an upper BVH query
  /// limit, then asks the facet BVH for a squared safety estimate.
  /// @tparam Real_v Floating-point scalar type.
  /// @param tessellated Runtime tessellated data.
  /// @param point Local point.
  /// @param[out] safety Non-negative exit safety.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SafetyToOutNoConv(UnplacedStruct_t const &tessellated,
                                                                             Vector3D<Real_v> const &point,
                                                                             Real_v &safety)
  {
    int isurf;
    // get a quick upper limit from the min-distance to fixed set of anchor points on the surface
    // --> this limits BVH search from the start
    // TODO: this should also be vectorizable on the CPU

    float upper_limit_sq = InfinityLength<float>();
    const float px = point.x(), py = point.y(), pz = point.z();
    for (int i = 0; i < TessellatedRuntimeStruct<float>::N; ++i) {
      float dx       = tessellated.fTestPoints_x[i] - px;
      float dy       = tessellated.fTestPoints_y[i] - py;
      float dz       = tessellated.fTestPoints_z[i] - pz;
      float dist2    = dx * dx + dy * dy + dz * dz;
      upper_limit_sq = Min(dist2, upper_limit_sq);
    }

    constexpr bool approxSafety = true;
    Real_v safetysq             = SafetySq<Real_v, approxSafety, double>(tessellated, point, isurf, upper_limit_sq);
    safety                      = Sqrt(safetysq);

    // Keep the NoConv helper unsigned; public wrappers own wrong-side conventions.
    // safety on boundary should be zero --> see ShapeTester
    if (safety < kTolerance) {
      safety = 0.;
    }
  }

  /// @brief Compute safety from an inside or surface point to leave.
  /// @details Performs public wrong-side classification before delegating to
  /// the no-convention helper. Outside starts return `-1`, and surface starts
  /// return zero.
  /// @tparam Real_v Floating-point scalar type.
  /// @param tessellated Runtime tessellated data.
  /// @param point Local point.
  /// @param[out] safety Exit safety, `-1`, or zero for surface starts.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SafetyToOut(UnplacedStruct_t const &tessellated,
                                                                       Vector3D<Real_v> const &point, Real_v &safety)
  {
#ifndef VECGEOM_NO_WRONGSIDE_CONV
    Inside_t inside;
    Inside<Real_v>(tessellated, point, inside);
    if (inside == kOutside) {
      safety = Real_v(-1.);
      return;
    }
    // Surface-point safety is defined as zero; avoid the expensive BVH safety query in this case.
    if (inside == kSurface) {
      safety = Real_v(0.);
      return;
    }
#endif
    SafetyToOutNoConv<Real_v>(tessellated, point, safety);
  }

  /// @brief Compute a surface-facet normal, falling back to the closest facet.
  /// @details Accept a contained plane projection within the surface band;
  /// unresolved edge bands and off-surface points use the closest facet.
  /// @tparam Real_v Floating-point scalar type.
  /// @param tessellated Runtime tessellated data.
  /// @param point Local point.
  /// @param[out] valid True when a facet normal was identified.
  /// @return Facet normal, or zero when @p valid is false.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static Vector3D<Real_v> NormalKernel(UnplacedStruct_t const &tessellated,
                                                                                    Vector3D<Real_v> const &point,
                                                                                    bool &valid)
  {
    // Computes the normal on a surface and returns it as a unit vector
    int isurf = -1;

    // A contained plane projection in the surface band fixes a valid normal.
    // Query-box overlap prunes nodes; unresolved edge bands retain the nearest
    // search below without burdening the ordinary face test.
    // Round outward so narrowing cannot shrink the query surface band.
    Vector3D<float> lower(point - Vector3D<Real_v>(kHalfTolerance));
    Vector3D<float> upper(point + Vector3D<Real_v>(kHalfTolerance));
    for (int axis = 0; axis < 3; ++axis) {
      lower[axis] = std::nextafter(lower[axis], -InfinityLength<float>());
      upper[axis] = std::nextafter(upper[axis], InfinityLength<float>());
    }
    const auto *nodes  = tessellated.fBVH->GetNodes();
    const auto *primId = tessellated.fBVH->GetPrimId();
    unsigned int stack[FacetBVH_t::kTraversalStackSize];
    const auto overlaps = [&](unsigned int index) {
      const auto &bounds = nodes[index].bounds;
      // Combine the six slab tests as integers: branchless, and free of the
      // clang -Wbitwise-instead-of-logical warning raised by '&' on booleans.
      return (static_cast<int>(upper[0] >= bounds[0]) & static_cast<int>(lower[0] <= bounds[1]) &
              static_cast<int>(upper[1] >= bounds[2]) & static_cast<int>(lower[1] <= bounds[3]) &
              static_cast<int>(upper[2] >= bounds[4]) & static_cast<int>(lower[2] <= bounds[5])) != 0;
    };
    int size = 0;
    if (overlaps(0)) stack[size++] = 0;
    while (size) {
      const auto &node = nodes[stack[--size]];
      if (node.IsLeaf()) {
        for (int i = 0; i < node.count; ++i) {
          const auto &facet = tessellated.fFacets[primId[node.first + i]];
          if (IsOnFiniteFacet<Real_v, false>(facet, point)) {
            valid = true;
            return facet.fNormal;
          }
        }
      } else {
        // Reject children before stacking to avoid divergent visits to misses.
        if (overlaps(node.first + 1)) stack[size++] = node.first + 1;
        if (overlaps(node.first)) stack[size++] = node.first;
      }
    }

    SafetySq<Real_v, false, double>(tessellated, point, isurf);
    if (isurf != -1) {
      valid = true;
      return tessellated.fFacets[isurf].fNormal;
    }
    valid = false;
    return Vector3D<Real_v>(0., 0., 0.);
  }

  /// @brief Query the BVH for squared distance to the closest relevant facet.
  /// @tparam Real_v Floating-point scalar type returned to callers.
  /// @tparam ToIn Selects the BVH side convention for entry or exit safety.
  /// @tparam T Internal precision used for facet safety comparisons.
  /// @param tessellated Runtime tessellated data.
  /// @param point Local point.
  /// @param[out] isurf Index of the closest visited facet, or `-1`.
  /// @param limit_sq Initial squared-distance limit for BVH pruning.
  /// @return Squared safety estimate returned by the BVH query.
  template <typename Real_v, bool ToIn, typename T = float>
  VECCORE_ATT_HOST_DEVICE static Real_v SafetySq(UnplacedStruct_t const &tessellated, Vector3D<Real_v> const &point,
                                                 int &isurf, Real_v limit_sq = InfinityLength<Real_v>())
  {
    T safetysq = limit_sq;
    isurf      = -1;
    Vector3D<T> pointv(point);

    const auto safety_estimated_sq =
        VisitPointSafetyBoxes<ToIn>(*tessellated.fBVH, point, limit_sq, [&](int primID, float &bvhSafetySqr) {
          const T this_safety_sq = tessellated.fFacets[primID].template SafetySq<T>(pointv);
          if (this_safety_sq < safetysq) {
            safetysq = this_safety_sq;
            isurf    = primID;
          }
          bvhSafetySqr = vecCore::math::Min(bvhSafetySqr, static_cast<float>(this_safety_sq));
        });
    // if (isurf != -1) {
    //  this is the best value
    //   return safetysq;
    // }
    // Somehow we'll need to know the safety was estimated or correct

    return safety_estimated_sq; // an estimate
  }
}; // end TessellatedImplementation

} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif // VECGEOM_VOLUMES_KERNEL_TESSELLATEDIMPLEMENTATION_H_
