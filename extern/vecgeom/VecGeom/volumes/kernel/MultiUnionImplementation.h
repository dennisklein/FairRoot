//===-- kernel/MultiUnionImplementation.h ---------------------------*- C++ -*-===//
//===--------------------------------------------------------------------------===//
/// @file MultiUnionImplementation.h
/// @brief Navigation kernels for an IndexedBVH-backed union of placed components.
/// @author Mihaela Gheata (mihaela.gheata@cern.ch)

#ifndef VECGEOM_VOLUMES_KERNEL_MULTIUNIONIMPLEMENTATION_H_
#define VECGEOM_VOLUMES_KERNEL_MULTIUNIONIMPLEMENTATION_H_

#include "VecGeom/base/Vector3D.h"
#include "VecGeom/volumes/MultiUnionStruct.h"
#include "VecGeom/volumes/SurfaceHitView.h"
#include "VecGeom/volumes/kernel/GenericKernels.h"
#include <VecCore/VecCore>

#include <cstdio>

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE(struct MultiUnionImplementation;);
VECGEOM_DEVICE_DECLARE_CONV(struct, MultiUnionImplementation);

inline namespace VECGEOM_IMPL_NAMESPACE {

class PlacedMultiUnion;
struct MultiUnionStruct;
class UnplacedMultiUnion;

struct MultiUnionImplementation {

  using PlacedShape_t    = PlacedMultiUnion;
  using UnplacedStruct_t = MultiUnionStruct;
  using UnplacedVolume_t = UnplacedMultiUnion;

  /**
   * Visit component boxes containing @p point until the visitor accepts one.
   * @param munion MultiUnion runtime data.
   * @param point Query point in MultiUnion coordinates.
   * @param visit Exact component predicate returning true to stop traversal.
   * @return Whether the visitor accepted a component.
   */
  template <typename Real_v, typename Visitor>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static bool VisitPointCandidates(UnplacedStruct_t const &munion,
                                                                                Vector3D<Real_v> const &point,
                                                                                Visitor &&visit)
  {
#if defined(__clang__) && !defined(VECCORE_CUDA_DEVICE_COMPILATION)
// Prevent Clang from widening three float operations into a fourth poison lane under trapped FPEs.
#pragma STDC FENV_ACCESS ON
#pragma clang fp exceptions(strict)
#endif
    using BVH_t        = UnplacedStruct_t::BVHStructure;
    using BVHReal_t    = UnplacedStruct_t::BVHReal_t;
    auto const &bvh    = *munion.fNavHelper;
    auto const *nodes  = bvh.GetNodes();
    auto const *primId = bvh.GetPrimId();
    const Vector3D<BVHReal_t> bvhPoint(static_cast<BVHReal_t>(point[0]), static_cast<BVHReal_t>(point[1]),
                                       static_cast<BVHReal_t>(point[2]));
    if (!bvh.NodeContains(0, bvhPoint)) return false;

    unsigned int stack[BVH_t::kTraversalStackSize];
    int stackSize    = 0;
    unsigned int top = 0;
    for (;;) {
      while (!nodes[top].IsLeaf()) {
        const unsigned int childL = nodes[top].first;
        const unsigned int childR = childL + 1;
        const bool containsL      = bvh.NodeContains(childL, bvhPoint);
        const bool containsR      = bvh.NodeContains(childR, bvhPoint);
        if (containsL) {
          if (containsR) stack[stackSize++] = childR;
          top = childL;
        } else if (containsR) {
          top = childR;
        } else {
          goto pop;
        }
      }

      for (int i = 0; i < nodes[top].count; ++i) {
        const int boxId = primId[nodes[top].first + i];
        if (bvh.BoxContains(boxId, bvhPoint) && visit(static_cast<size_t>(bvh.GetDaughterId(boxId)))) return true;
      }

    pop:
      if (stackSize == 0) return false;
      top = stack[--stackSize];
    }
  }

  /**
   * Traverse component boxes near-first while the visitor reduces the ray bound.
   * @param munion MultiUnion runtime data.
   * @param point Ray origin in MultiUnion coordinates.
   * @param direction Normalized ray direction.
   * @param stepMax Maximum query distance.
   * @param visit Exact component query that may reduce the BVH-precision distance bound.
   */
  template <typename Real_v, typename Visitor>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void VisitRayCandidates(UnplacedStruct_t const &munion,
                                                                              Vector3D<Real_v> const &point,
                                                                              Vector3D<Real_v> const &direction,
                                                                              Real_v const &stepMax, Visitor &&visit)
  {
#if defined(__clang__) && !defined(VECCORE_CUDA_DEVICE_COMPILATION)
// Prevent Clang from widening three float operations into a fourth poison lane under trapped FPEs.
#pragma STDC FENV_ACCESS ON
#pragma clang fp exceptions(strict)
#endif
    using BVH_t        = UnplacedStruct_t::BVHStructure;
    using BVHReal_t    = UnplacedStruct_t::BVHReal_t;
    auto const &bvh    = *munion.fNavHelper;
    auto const *nodes  = bvh.GetNodes();
    auto const *primId = bvh.GetPrimId();
    const Vector3D<BVHReal_t> bvhPoint(static_cast<BVHReal_t>(point[0]), static_cast<BVHReal_t>(point[1]),
                                       static_cast<BVHReal_t>(point[2]));
    const Vector3D<BVHReal_t> invDirection(BVHReal_t(1.) / vecgeom::NonZero(static_cast<BVHReal_t>(direction[0])),
                                           BVHReal_t(1.) / vecgeom::NonZero(static_cast<BVHReal_t>(direction[1])),
                                           BVHReal_t(1.) / vecgeom::NonZero(static_cast<BVHReal_t>(direction[2])));
    const unsigned int octant = (invDirection[0] < BVHReal_t(0.) ? 1u : 0u) |
                                (invDirection[1] < BVHReal_t(0.) ? 2u : 0u) |
                                (invDirection[2] < BVHReal_t(0.) ? 4u : 0u);
    BVHReal_t bvhStep         = static_cast<BVHReal_t>(stepMax);
    BVHReal_t approach;
    if (!bvh.IntersectNode(0, bvhPoint, invDirection, octant, bvhStep, approach)) return;

    unsigned int stack[BVH_t::kTraversalStackSize];
    int stackSize    = 0;
    unsigned int top = 0;
    for (;;) {
      while (!nodes[top].IsLeaf()) {
        const unsigned int childL = nodes[top].first;
        const unsigned int childR = childL + 1;
        BVHReal_t approachL, approachR;
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
        const int boxId = primId[nodes[top].first + i];
        if (!bvh.IntersectBox(boxId, bvhPoint, invDirection, bvhStep, approach)) continue;
        visit(static_cast<size_t>(bvh.GetDaughterId(boxId)), bvhStep);
      }

    pop:
      if (stackSize == 0) return;
      top = stack[--stackSize];
    }
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void Contains(UnplacedStruct_t const &munion,
                                                                    Vector3D<Real_v> const &point, bool &inside)
  {
    inside = false;
    VisitPointCandidates(munion, point, [&](size_t id) {
      inside = munion.fVolumes[id]->Contains(point);
      return inside;
    });
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void Inside(UnplacedStruct_t const &munion,
                                                                  Vector3D<Real_v> const &point, Inside_t &inside)
  {
    inside                = EInside::kOutside;
    unsigned surfaceCount = 0;
    VisitPointCandidates(munion, point, [&](size_t id) {
      auto inside_crt = munion.fVolumes[id]->Inside(point);
      if (inside_crt == EInside::kInside) {
        inside = EInside::kInside;
        return true;
      }
      if (inside_crt == EInside::kSurface) {
        inside = EInside::kSurface;
        ++surfaceCount;
      }

      return false;
    });
    if (inside != EInside::kSurface || surfaceCount < 2) return;

    // Match Boolean-union classification at internal contacts. Only multiple
    // surface candidates need normals; ordinary inside/outside queries stop above.
    VisitPointCandidates(munion, point, [&](size_t id) {
      auto const *volume = munion.fVolumes[id];
      if (volume->Inside(point) != EInside::kSurface) return false;
      Vector3D<Real_v> normal;
      volume->Normal(point, normal);
      const bool internal = VisitPointCandidates(munion, point, [&](size_t other) {
        if (other >= id || munion.fVolumes[other]->Inside(point) != EInside::kSurface) return false;
        Vector3D<Real_v> otherNormal;
        munion.fVolumes[other]->Normal(point, otherNormal);
        return normal.Dot(otherNormal) < 0;
      });
      if (internal) inside = EInside::kInside;
      return internal;
    });
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void InsideComponent(UnplacedStruct_t const &munion,
                                                                           Vector3D<Real_v> const &point,
                                                                           int &component)
  {
    component = -1;
    VisitPointCandidates(munion, point, [&](size_t id) {
      auto inside_crt = munion.fVolumes[id]->Inside(point);
      if (inside_crt != EInside::kOutside) {
        component = id;
        return true;
      }
      return false;
    });
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void InsideCluster(UnplacedStruct_t const &munion,
                                                                         Vector3D<Real_v> const &point,
                                                                         Vector3D<Real_v> const &direction,
                                                                         int &component)
  {
    // loop cluster overlap candidates for current component, then update component
    const size_t begin = munion.fNeighbourOffsets[component];
    const size_t end   = munion.fNeighbourOffsets[component + 1];
    for (size_t i = begin; i < end; ++i) {
      const size_t neighbour = munion.fNeighbours[i];
      const auto inside      = munion.fVolumes[neighbour]->Inside(point);
      if (inside == EInside::kOutside) continue;
      if (inside == EInside::kSurface) {
        // A shallow handoff can remain in the surface tolerance band after
        // the push. Accept inward/tangent motion, but do not relocate through
        // a shared exterior face when the ray is leaving both components.
        Vector3D<Real_v> normal;
        munion.fVolumes[neighbour]->Normal(point, normal);
        if (normal.Dot(direction) > 0) continue;
      }
      component = neighbour;
      return;
    }
    component = -1;
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void DistanceToIn(UnplacedStruct_t const &munion,
                                                                        Vector3D<Real_v> const &point,
                                                                        Vector3D<Real_v> const &direction,
                                                                        Real_v const &stepMax, Real_v &distance)
  {
    distance = kInfLength;
    VisitRayCandidates(
        munion, point, direction, stepMax, [&](size_t id, typename UnplacedStruct_t::BVHReal_t &bvhStep) {
          const Real_v distance_crt = munion.fVolumes[id]->DistanceToIn(point, direction, Min(stepMax, distance));
          if (distance_crt < distance) {
            distance = distance_crt;
            bvhStep  = static_cast<typename UnplacedStruct_t::BVHReal_t>(distance);
          }
        });
  }

  /**
   * @brief Publish a component plane only when it supports the complete union.
   * @param munion Component data and padded union bounds.
   * @param hitPoint Component exit point in union coordinates.
   * @param component Index of the exiting component.
   * @param componentHit Component supporting-plane result.
   * @param componentNormal Component-local normal, when requested.
   * @param[out] hitInfo Requested union exit information.
   * @details The union AABB gives a constant-cost sufficient test; rejection is conservative.
   */
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SetExitHit(UnplacedStruct_t const &munion,
                                                                      Vector3D<Real_v> const &hitPoint,
                                                                      size_t component,
                                                                      SurfaceHitView<Real_v> const &componentHit,
                                                                      Vector3D<Real_v> const &componentNormal,
                                                                      SurfaceHitView<Real_v> *hitInfo)
  {
    bool leftBehind = componentHit.fLeftBehind;
    Vector3D<Real_v> normal;
    if (leftBehind && (munion.fVolumes.size() > 1 || hitInfo->WantsNormal())) {
      normal = munion.fVolumes[component]->GetTransformation()->InverseTransformDirection(componentNormal);
    }
    if (leftBehind && munion.fVolumes.size() > 1) {
      // The component normal supports the complete union when its plane also
      // supports the union AABB. Account for the per-axis AABB padding; a
      // failed test is conservative.
      Vector3D<Real_v> const support(normal.x() >= 0 ? munion.fMaxExtent.x() : munion.fMinExtent.x(),
                                     normal.y() >= 0 ? munion.fMaxExtent.y() : munion.fMinExtent.y(),
                                     normal.z() >= 0 ? munion.fMaxExtent.z() : munion.fMinExtent.z());
      Real_v const padding = kAABBPadding * (Abs(normal.x()) + Abs(normal.y()) + Abs(normal.z()));
      leftBehind           = normal.Dot(support - hitPoint) <= padding + kTolerance;
    }
    hitInfo->SetLeftBehind(leftBehind);
    if (leftBehind && hitInfo->WantsNormal()) hitInfo->SetNormal(normal);
  }

  /// Execute the existing traversal with hit recording compiled in or out.
  template <typename Real_v, bool RecordHitInfo>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void DistanceToOutImpl(UnplacedStruct_t const &munion,
                                                                             Vector3D<Real_v> const &point,
                                                                             Vector3D<Real_v> const &direction,
                                                                             Real_v const &stepMax, Real_v &distance,
                                                                             SurfaceHitView<Real_v> *hitInfo)
  {
    constexpr Real_v eps = 10 * kTolerance;
    distance             = -1.;
    if constexpr (RecordHitInfo) hitInfo->Clear();
    // Locate the component containing the point
    int comp;
    InsideComponent(munion, point, comp);
    if (comp < 0) return; // Point not inside
    // Compute distance to exit current component
    distance              = -eps;
    Real_v dstep          = 1.;
    Vector3D<Real_v> pnew = point;
    Vector3D<Real_v> local, ldir;
    // A zero exit at an internal contact can hand off to another component.
    // Terminate when relocation finds no neighbour, not on that first zero step.
    while (comp >= 0) {
      size_t component = (size_t)comp;
      munion.fVolumes[component]->GetTransformation()->Transform(pnew, local);
      munion.fVolumes[component]->GetTransformation()->TransformDirection(direction, ldir);
      Vector3D<Real_v> componentNormal;
      SurfaceHitView<Real_v> componentHit;
      if constexpr (RecordHitInfo) {
        componentHit.fNormal = (munion.fVolumes.size() > 1 || hitInfo->WantsNormal()) ? &componentNormal : nullptr;
        dstep                = munion.fVolumes[component]->DistanceToOut(local, ldir, stepMax, &componentHit);
      } else {
        dstep = munion.fVolumes[component]->DistanceToOut(local, ldir, stepMax);
      }
      VECGEOM_ASSERT(dstep < kInfLength);
      distance += dstep + eps;
      // If no neighbours, exit
      if (munion.fNeighbourOffsets[component] == munion.fNeighbourOffsets[component + 1]) {
        if constexpr (RecordHitInfo)
          SetExitHit(munion, pnew + dstep * direction, component, componentHit, componentNormal, hitInfo);
        return;
      }
      // Propagate to exit of current component
      pnew += (dstep + eps) * direction;
      // Try to relocate inside the cluster of neighbours
      MultiUnionImplementation::InsideCluster(munion, pnew, direction, comp);
      if constexpr (RecordHitInfo) {
        if (comp < 0) SetExitHit(munion, pnew - eps * direction, component, componentHit, componentNormal, hitInfo);
      }
    }
  }

  // A component may itself be a MultiUnion. Keep a non-forced-inline boundary
  // on the hit-aware path so speculative devirtualization can stop expanding it.
  template <typename Real_v>
  VECCORE_ATT_HOST_DEVICE static void DistanceToOutWithHit(UnplacedStruct_t const &munion,
                                                           Vector3D<Real_v> const &point,
                                                           Vector3D<Real_v> const &direction, Real_v const &stepMax,
                                                           Real_v &distance, SurfaceHitView<Real_v> *hitInfo)
  {
    DistanceToOutImpl<Real_v, true>(munion, point, direction, stepMax, distance, hitInfo);
  }

  /**
   * @brief Return the union exit distance and optional supporting-plane normal.
   * @tparam Real_v Scalar floating-point type.
   * @param munion Component data and navigation hierarchy.
   * @param point Starting point in union coordinates.
   * @param direction Unit ray direction.
   * @param stepMax Distance bound forwarded with the existing component convention.
   * @param[out] distance Exit distance, or the existing wrong-side sentinel.
   * @param[out] hitInfo Optional output; false leaves caller normal storage untouched.
   * @details A component plane is published only if it also supports the union.
   * No hit recording or support test is performed when hitInfo is null.
   */
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void DistanceToOut(UnplacedStruct_t const &munion,
                                                                         Vector3D<Real_v> const &point,
                                                                         Vector3D<Real_v> const &direction,
                                                                         Real_v const &stepMax, Real_v &distance,
                                                                         SurfaceHitView<Real_v> *hitInfo = nullptr)
  {
    if (hitInfo)
      DistanceToOutWithHit(munion, point, direction, stepMax, distance, hitInfo);
    else
      DistanceToOutImpl<Real_v, false>(munion, point, direction, stepMax, distance, nullptr);
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SafetyToInComp(UnplacedStruct_t const &munion,
                                                                          Vector3D<Real_v> const &point, Real_v &safety,
                                                                          int &component)
  {
    using BVH_t        = UnplacedStruct_t::BVHStructure;
    using BVHReal_t    = UnplacedStruct_t::BVHReal_t;
    auto const &bvh    = *munion.fNavHelper;
    auto const *nodes  = bvh.GetNodes();
    auto const *primId = bvh.GetPrimId();
    const Vector3D<BVHReal_t> bvhPoint(static_cast<BVHReal_t>(point[0]), static_cast<BVHReal_t>(point[1]),
                                       static_cast<BVHReal_t>(point[2]));
    unsigned int stack[BVH_t::kTraversalStackSize];
    int stackSize = 1;
    stack[0]      = 0;
    safety        = vecgeom::InfinityLength<Real_v>();
    component     = -1;

    while (stackSize > 0) {
      const unsigned int top = stack[--stackSize];
      if (static_cast<Real_v>(bvh.NodeSafetySqr(top, bvhPoint)) > safety * safety) continue;
      if (nodes[top].IsLeaf()) {
        for (int i = 0; i < nodes[top].count; ++i) {
          const int boxId = primId[nodes[top].first + i];
          if (bvh.BoxSafety(boxId, point) >= safety) continue;
          const int id           = bvh.GetDaughterId(boxId);
          const Real_v safetycrt = munion.fVolumes[id]->SafetyToIn(point);
          if (safetycrt > 0 && safetycrt < safety) {
            safety    = safetycrt;
            component = id;
          }
        }
        continue;
      }

      const unsigned int childL = nodes[top].first;
      const unsigned int childR = childL + 1;
      const BVHReal_t safetyL   = bvh.NodeSafetySqr(childL, bvhPoint);
      const BVHReal_t safetyR   = bvh.NodeSafetySqr(childR, bvhPoint);
      const Real_v limitSqr     = safety * safety;
      const bool traverseL      = static_cast<Real_v>(safetyL) < limitSqr;
      const bool traverseR      = static_cast<Real_v>(safetyR) < limitSqr;
      if (safetyR < safetyL) {
        if (traverseL) stack[stackSize++] = childL;
        if (traverseR) stack[stackSize++] = childR;
      } else {
        if (traverseR) stack[stackSize++] = childR;
        if (traverseL) stack[stackSize++] = childL;
      }
    }
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SafetyToIn(UnplacedStruct_t const &munion,
                                                                      Vector3D<Real_v> const &point, Real_v &safety)
  {
    int comp;
    InsideComponent(munion, point, comp);
    if (comp > -1) {
      safety = -1.;
      return;
    }

    SafetyToInComp<Real_v>(munion, point, safety, comp);
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SafetyToOut(UnplacedStruct_t const &munion,
                                                                       Vector3D<Real_v> const &point, Real_v &safety)
  {
    // Locate the component containing the point
    int comp;
    MultiUnionImplementation::InsideComponent(munion, point, comp);
    if (comp < 0) {
      safety = -1.; // Point not inside
      return;
    }
    // Compute safety to exit current component
    Vector3D<Real_v> const local = munion.fVolumes[comp]->GetTransformation()->Transform(point);
    safety                       = munion.fVolumes[comp]->SafetyToOut(local);
    VECGEOM_ASSERT(safety > -kTolerance);
    // Loop cluster of neighbours
    const size_t begin = munion.fNeighbourOffsets[comp];
    const size_t end   = munion.fNeighbourOffsets[comp + 1];
    for (size_t i = begin; i < end; ++i) {
      const size_t neighbour       = munion.fNeighbours[i];
      Vector3D<Real_v> const local = munion.fVolumes[neighbour]->GetTransformation()->Transform(point);
      Real_v safetycrt             = munion.fVolumes[neighbour]->SafetyToOut(local);
      if (safetycrt > 0 && safetycrt < safety) safety = safetycrt;
    }
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static Vector3D<Real_v> NormalKernel(UnplacedStruct_t const &munion,
                                                                                    Vector3D<Real_v> const &point,
                                                                                    bool &valid)
  {
    // Locate the component containing the point
    int comp;
    valid = false;
    Vector3D<Real_v> direction;

    InsideComponent(munion, point, comp);
    // If component not found, locate closest one
    Real_v safety;
    if (comp < 0) SafetyToInComp(munion, point, safety, comp);
    if (comp < 0) return direction;

    Vector3D<Real_v> local = munion.fVolumes[comp]->GetTransformation()->Transform(point);
    Vector3D<Real_v> ldir;
    valid = munion.fVolumes[comp]->GetUnplacedVolume()->Normal(local, ldir);
    if (valid) direction = munion.fVolumes[comp]->GetTransformation()->InverseTransformDirection(ldir);
    return direction;
  }

}; // end MultiUnionImplementation
} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif // VECGEOM_VOLUMES_KERNEL_MULTIUNIONIMPLEMENTATION_H_
