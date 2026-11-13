// SPDX-FileCopyrightText: 2020 CERN
// SPDX-License-Identifier: Apache-2.0

/**
 * @file BVHNavigator.h
 * @brief Navigation methods for geometry.
 */

#ifndef BVH_NAVIGATOR_H
#define BVH_NAVIGATOR_H

#include "VecGeom/base/Global.h"
#include "VecGeom/base/Math.h"
#include "VecGeom/base/Vector3D.h"
#include "VecGeom/navigation/NavigationState.h"
#include "VecGeom/volumes/LogicalVolume.h"
#include "VecGeom/management/BVHManager.h"
#include "VecGeom/management/GeoManager.h"

#ifdef VECGEOM_ENABLE_CUDA
#include "VecGeom/backend/cuda/Interface.h"
#endif

#include <type_traits>

namespace vgbrep {
namespace protonav {
template <typename Real_t>
class BVHSurfNavigator;
}
} // namespace vgbrep

namespace vecgeom {

inline namespace VECGEOM_IMPL_NAMESPACE {
class BVHSafetyEstimator;
}

/**
 * @brief IndexedBVH-accelerated navigation for placed-volume hierarchies.
 *
 * Application code uses the point-location and tracking interfaces documented
 * below. In particular, AdePT uses `ComputeSafety()`,
 * `ComputeStepAndNextVolume()`, and `RelocateToNextVolume()`. The component and
 * diagnostic interfaces are intended for VecGeom navigation implementations
 * and tests, not for transport applications.
 */
class BVHNavigator {

public:
  /** @name Application navigation interface
   *  Transport clients use these methods; BVH lookup and traversal remain internal.
   *  @{ */

  /// Boundary displacement used before or after navigation across a surface.
  static constexpr Precision kBoundaryPush = 10 * vecgeom::kTolerance;

  /**
   * @brief Locate the deepest daughter containing a point.
   * @param vol Starting placed volume.
   * @param point Point in the parent frame of @p vol.
   * @param[in,out] path State initially at the parent of @p vol; extended to the deepest containing daughter.
   * @param top Whether to first require that @p vol contains @p point.
   * @param exclude Optional volume excluded from the first daughter query.
   * @return Deepest containing volume, or null when the optional top check fails.
   */
  VECCORE_ATT_HOST_DEVICE
  static VECGEOM_FORCE_INLINE Daughter LocatePointIn(vecgeom::VPlacedVolume const *vol,
                                                     Vector3D<Precision> const &point, vecgeom::NavigationState &path,
                                                     bool top, vecgeom::VPlacedVolume const *exclude = nullptr)
  {
    return LocatePointInImpl(vol, point, path, top, exclude);
  }

  /**
   * @brief Relocate a local point after leaving the current top volume.
   * @param localpoint Point in the local frame of `path.Top()`.
   * @param[in,out] path State to climb and descend to the new containing volume.
   * @return Deepest containing volume, or null after leaving the world.
   */
  VECCORE_ATT_HOST_DEVICE
  static VECGEOM_FORCE_INLINE Daughter RelocatePoint(Vector3D<Precision> const &localpoint,
                                                     vecgeom::NavigationState &path)
  {
    return RelocatePointImpl(localpoint, path);
  }

  /**
   * @brief Locate the deepest daughter from the current navigation state.
   * @param localpoint Point in the local frame of `path.Top()`.
   * @param[in,out] path Starting state, extended to the deepest containing daughter.
   * @param top Whether to first require that the current top contains @p localpoint.
   * @param exclude Optional volume excluded from the first daughter query.
   * @return Final top volume, or null when the optional top check fails.
   */
  VECCORE_ATT_HOST_DEVICE
  static VECGEOM_FORCE_INLINE Daughter LocatePointInNavState(Vector3D<Precision> const &localpoint,
                                                             vecgeom::NavigationState &path, bool top,
                                                             vecgeom::VPlacedVolume const *exclude = nullptr)
  {
    return LocatePointInNavStateImpl(localpoint, path, top, exclude);
  }

  /**
   * @brief Compute isotropic Safety from a global point in the current state.
   * @param globalpoint Point in global coordinates.
   * @param state Current navigation state.
   * @param limit Accuracy limit; values above it may be conservative.
   * @return Safety to the nearest boundary, or infinity when @p state is outside.
   */
  VECCORE_ATT_HOST_DEVICE
  static VECGEOM_FORCE_INLINE Precision ComputeSafety(Vector3D<Precision> const &globalpoint,
                                                      vecgeom::NavigationState const &state,
                                                      Precision limit = InfinityLength<Precision>())
  {
    return ComputeSafetyImpl(globalpoint, state, limit);
  }

  /**
   * @brief Compute a bounded step and fully relocate the output state.
   * @param globalpoint Point in global coordinates and inside @p in_state.
   * @param globaldir Unit direction in global coordinates.
   * @param step_limit Maximum requested step.
   * @param in_state Current navigation state.
   * @param[out] out_state State at the endpoint, including transitive boundary relocation.
   * @return Geometry-limited or requested step length.
   */
  VECCORE_ATT_HOST_DEVICE
  static VECGEOM_FORCE_INLINE Precision ComputeStepAndPropagatedState(Vector3D<Precision> const &globalpoint,
                                                                      Vector3D<Precision> const &globaldir,
                                                                      Precision step_limit,
                                                                      vecgeom::NavigationState const &in_state,
                                                                      vecgeom::NavigationState &out_state)
  {
    return ComputeStepAndPropagatedStateImpl(globalpoint, globaldir, step_limit, in_state, out_state);
  }

  /**
   * @brief Compute a bounded step and update the state to the first crossed boundary.
   *
   * The output state enters the hit daughter or leaves the current volume, but
   * does not resolve additional coincident boundaries. Call
   * `RelocateToNextVolume()` after moving the track to the returned endpoint.
   *
   * @param globalpoint Point in global coordinates and inside @p in_state.
   * @param globaldir Unit direction in global coordinates.
   * @param step_limit Maximum requested step.
   * @param in_state Current navigation state.
   * @param[out] out_state State at the first crossed boundary, or a copy of @p in_state for a physics step.
   * @return Geometry-limited or requested step length.
   */
  VECCORE_ATT_HOST_DEVICE
  static VECGEOM_FORCE_INLINE Precision ComputeStepAndNextVolume(Vector3D<Precision> const &globalpoint,
                                                                 Vector3D<Precision> const &globaldir,
                                                                 Precision step_limit,
                                                                 vecgeom::NavigationState const &in_state,
                                                                 vecgeom::NavigationState &out_state)
  {
    return ComputeStepAndNextVolumeImpl(globalpoint, globaldir, step_limit, in_state, out_state);
  }

  /**
   * @brief Compute a bounded approach step without changing the navigation state.
   * @param globalpoint Point in global coordinates and inside @p in_state.
   * @param globaldir Unit direction in global coordinates.
   * @param step_limit Maximum requested step.
   * @param in_state Current navigation state.
   * @return Step limited by the next daughter approach or by @p step_limit.
   */
  VECCORE_ATT_HOST_DEVICE
  static VECGEOM_FORCE_INLINE Precision ComputeStepToApproachNextVolume(Vector3D<Precision> const &globalpoint,
                                                                        Vector3D<Precision> const &globaldir,
                                                                        Precision step_limit,
                                                                        vecgeom::NavigationState const &in_state)
  {
    return ComputeStepToApproachNextVolumeImpl(globalpoint, globaldir, step_limit, in_state);
  }

  /**
   * @brief Complete relocation after `ComputeStepAndNextVolume()`.
   * @param globalpoint Track endpoint in global coordinates.
   * @param globaldir Track direction in global coordinates.
   * @param[in,out] state Boundary state to descend through coincident volumes.
   */
  VECCORE_ATT_HOST_DEVICE
  static VECGEOM_FORCE_INLINE void RelocateToNextVolume(Vector3D<Precision> const &globalpoint,
                                                        Vector3D<Precision> const &globaldir,
                                                        vecgeom::NavigationState &state)
  {
    RelocateToNextVolumeImpl(globalpoint, globaldir, state);
  }

  /** @} */

private:
  template <typename Real_t>
  friend class ::vgbrep::protonav::BVHSurfNavigator;
  friend class BVHSafetyEstimator;

// to avoid changing the user code via templates, the precision is chosen at compile-time in the BVHManager and
// BVHNavigator
#ifdef VECGEOM_BVH_SINGLE
  using Real_t = float;
#else
  using Real_t = double;
#endif

  using Self_t = BVHNavigator;

  /// Host Safety scans up to this many boxes directly to avoid stack and internal-node overhead.
  static constexpr int kDirectSafetyBoxThreshold = 8;

  /** Box-traversal counters collected by the instrumented daughter ray query. */
  struct IntersectStatistics {
    long primitiveBoxTests{0};     ///< Primitive AABBs tested in visited leaves.
    long visitedLeaves{0};         ///< Leaf nodes reached by traversal.
    long cutNodes{0};              ///< Root-node rejections.
    long scheduledNodes{0};        ///< Child nodes whose AABBs intersect the ray.
    long intersectedBoxes{0};      ///< Primitive AABBs intersecting the ray.
    long deferredNodes{0};         ///< Far children deferred on the traversal stack.
    long prunableDeferredNodes{0}; ///< Deferred children beyond a later reduced step.
    long maximumStackDepth{0};     ///< Maximum number of deferred nodes.
  };

  /** Box-traversal counters collected by the instrumented daughter Safety query. */
  struct SafetyStatistics {
    long visitedNodes{0};      ///< Nodes removed from the traversal stack.
    long visitedLeaves{0};     ///< Leaf nodes reached by traversal.
    long primitiveBoxTests{0}; ///< Primitive AABB Safety evaluations.
    long relevantBoxes{0};     ///< Primitive AABBs closer than the current bound.
    long cutNodes{0};          ///< Nodes rejected by the current Safety bound.
  };

  VECCORE_ATT_HOST_DEVICE
  static auto GetBVH(int logicalVolumeId) { return vecgeom::BVHManager::GetBVH(logicalVolumeId); }

  /// @tparam conditional_exclusion If true, a daughter at a non-positive distance is refused when it is the excluded
  /// one or when the ray is leaving it.
  template <bool collect_statistics, typename Navigator, bool conditional_exclusion = false, typename BVHReal_t,
            typename Real_i>
  VECCORE_ATT_HOST_DEVICE static void CheckDaughterIntersectionsImpl(
      const vecgeom::IndexedBVH<BVHReal_t> &bvh, const Vector3D<Real_i> &localpoint, const Vector3D<Real_i> &localdir,
      Real_i &step, long const last_exited_id, long &hitcandidate_index, IntersectStatistics *statistics = nullptr,
      long *candidate_calls = nullptr)
  {
    struct StackItem {
      unsigned int node;
      BVHReal_t approach;
    };
    using StackEntry = std::conditional_t<collect_statistics, StackItem, unsigned int>;
    StackEntry stack[vecgeom::IndexedBVH<BVHReal_t>::kTraversalStackSize];
    int stackSize = 0;

    if constexpr (collect_statistics) {
      VECGEOM_ASSERT(statistics != nullptr);
      *statistics = IntersectStatistics{};
    }

    const bool mapped = bvh.HasDaughterMap();
    const Vector3D<BVHReal_t> point(static_cast<BVHReal_t>(localpoint[0]), static_cast<BVHReal_t>(localpoint[1]),
                                    static_cast<BVHReal_t>(localpoint[2]));
    const Vector3D<BVHReal_t> invDirection(static_cast<BVHReal_t>(1.) / vecgeom::NonZero(localdir[0]),
                                           static_cast<BVHReal_t>(1.) / vecgeom::NonZero(localdir[1]),
                                           static_cast<BVHReal_t>(1.) / vecgeom::NonZero(localdir[2]));
    // CUDA uses the direction-sign bits to select near/far slab endpoints.
    // The host predicate uses fixed endpoint loads; inlining removes this unused value there.
    const unsigned int octant = (invDirection[0] < BVHReal_t(0.) ? 1u : 0u) |
                                (invDirection[1] < BVHReal_t(0.) ? 2u : 0u) |
                                (invDirection[2] < BVHReal_t(0.) ? 4u : 0u);
    BVHReal_t bvhStep         = static_cast<BVHReal_t>(step);
    BVHReal_t rootApproach;
    if (!bvh.IntersectNode(0, point, invDirection, octant, bvhStep, rootApproach)) {
      if constexpr (collect_statistics) ++statistics->cutNodes;
      return;
    }

    auto const *nodes  = bvh.GetNodes();
    auto const *primId = bvh.GetPrimId();
    unsigned int top   = 0;
    for (;;) {
      while (!nodes[top].IsLeaf()) {
        const unsigned int childL = nodes[top].first;
        const unsigned int childR = childL + 1;
        BVHReal_t approachL, approachR;
        const bool hitL = bvh.IntersectNode(childL, point, invDirection, octant, bvhStep, approachL);
        const bool hitR = bvh.IntersectNode(childR, point, invDirection, octant, bvhStep, approachR);
        if constexpr (collect_statistics)
          statistics->scheduledNodes += static_cast<long>(hitL) + static_cast<long>(hitR);

        if (hitL) {
          // Visit the closer child first so exact hits reduce bvhStep before the deferred subtree.
          unsigned int near = childL;
          if (hitR) {
            unsigned int far   = childR;
            BVHReal_t farEntry = approachR;
            if (approachR < approachL) {
              near     = childR;
              far      = childL;
              farEntry = approachL;
            }
            if constexpr (collect_statistics) {
              stack[stackSize++] = StackItem{far, farEntry};
              ++statistics->deferredNodes;
              statistics->maximumStackDepth = vecgeom::Max(statistics->maximumStackDepth, long(stackSize));
            } else {
              stack[stackSize++] = far;
            }
          }
          top = near;
        } else if (hitR) {
          top = childR;
        } else {
          goto pop;
        }
      }

      if constexpr (collect_statistics) ++statistics->visitedLeaves;
      for (int i = 0; i < nodes[top].count; ++i) {
        if constexpr (collect_statistics) ++statistics->primitiveBoxTests;
        const int boxId = primId[nodes[top].first + i];
        BVHReal_t boxApproach;
        if (!bvh.IntersectBox(boxId, point, invDirection, bvhStep, boxApproach)) continue;
        if constexpr (collect_statistics) ++statistics->intersectedBoxes;

        const int prim      = bvh.GetDaughterId(boxId);
        const bool excluded = last_exited_id >= 0 && Navigator::SkipItem(bvh.GetRootId(), prim, last_exited_id);
        if constexpr (!conditional_exclusion) {
          if (excluded) continue;
        }
        if constexpr (collect_statistics) {
          if (candidate_calls) ++*candidate_calls;
        }

        const auto approach = static_cast<Real_i>(boxApproach);
        Real_i dist;
        if (!mapped) {
          dist =
              Navigator::CandidateDistanceToIn(bvh.GetRootId(), prim, localpoint + approach * localdir, localdir, step);
          if (dist > Real_i(0.)) dist += approach;
        } else if (conditional_exclusion) {
          // Take the distance from the original point: a later split box can start inside a refused daughter, or
          // in its hole, and the distance from there would accept it again.
          dist = Navigator::CandidateDistanceToIn(bvh.GetRootId(), prim, localpoint, localdir, step);
        } else {
          dist = Navigator::CandidateDistanceToIn(bvh.GetRootId(), prim, localpoint + approach * localdir, localdir,
                                                  vecgeom::Max(step - approach, Real_i(0.)));
          if (dist < Real_i(0.)) dist = Real_i(0.);
          dist += approach;
        }
        if (dist < step) {
          if constexpr (conditional_exclusion) {
            // As SimpleABBoxNavigator: at a non-positive distance refuse the excluded daughter and any daughter
            // the ray is leaving.
            if (dist <= Real_i(0.) &&
                (excluded || Navigator::CandidateLeaves(bvh.GetRootId(), prim, localpoint, localdir)))
              continue;
          }
          dist               = vecgeom::Max(dist, Real_i(0.));
          step               = dist;
          bvhStep            = static_cast<BVHReal_t>(dist);
          hitcandidate_index = prim;
        }
      }

    pop:
      if (stackSize == 0) return;
      if constexpr (collect_statistics) {
        const auto deferred = stack[--stackSize];
        if (deferred.approach > bvhStep) ++statistics->prunableDeferredNodes;
        top = deferred.node;
      } else {
        top = stack[--stackSize];
      }
    }
  }

  template <bool collect_statistics, typename Navigator, typename BVHReal_t>
  VECCORE_ATT_HOST_DEVICE static Precision ComputeBVHSafetyImpl(const vecgeom::IndexedBVH<BVHReal_t> &bvh,
                                                                Vector3D<Precision> const &localpoint, Precision safety,
                                                                Precision limit, Precision min_valid_safety,
                                                                SafetyStatistics *statistics = nullptr,
                                                                long *limit_shortcuts        = nullptr,
                                                                long *candidate_calls        = nullptr)
  {
    const Vector3D<BVHReal_t> point(static_cast<BVHReal_t>(localpoint[0]), static_cast<BVHReal_t>(localpoint[1]),
                                    static_cast<BVHReal_t>(localpoint[2]));
    if constexpr (collect_statistics) {
      VECGEOM_ASSERT(statistics != nullptr);
      *statistics = SafetyStatistics{};
    }

#ifndef VECCORE_CUDA_DEVICE_COMPILATION
    if constexpr (!collect_statistics) {
      // For small CPU hierarchies, one root rejection plus a direct box scan is
      // cheaper than internal-node tests and stack traffic. CUDA retains the tree path.
      if (bvh.GetRootNChild() <= kDirectSafetyBoxThreshold) {
        if (static_cast<Precision>(bvh.NodeSafetySqr(0, point)) > safety * safety) return safety;
        for (int index = 0; index < bvh.GetRootNChild(); ++index) {
          const int boxId           = bvh.GetPrimId()[index];
          const Precision boxSafety = bvh.BoxSafety(boxId, localpoint);
          if (boxSafety >= safety) continue;
          if (boxSafety > limit) {
            safety = vecgeom::Min(safety, boxSafety);
            continue;
          }
          const int prim       = bvh.GetDaughterId(boxId);
          const Precision dist = Navigator::CandidateSafetyToIn(bvh.GetRootId(), prim, localpoint);
          if (dist > min_valid_safety) safety = vecgeom::Min(safety, dist);
        }
        return safety;
      }
    }
#endif

    unsigned int stack[vecgeom::IndexedBVH<BVHReal_t>::kTraversalStackSize];
    int stackSize      = 1;
    stack[0]           = 0;
    auto const *nodes  = bvh.GetNodes();
    auto const *primId = bvh.GetPrimId();

    while (stackSize > 0) {
      const unsigned int id = stack[--stackSize];
      if constexpr (collect_statistics) ++statistics->visitedNodes;
      if (static_cast<Precision>(bvh.NodeSafetySqr(id, point)) > safety * safety) {
        if constexpr (collect_statistics) ++statistics->cutNodes;
        continue;
      }

      if (nodes[id].IsLeaf()) {
        if constexpr (collect_statistics) ++statistics->visitedLeaves;
        for (int i = 0; i < nodes[id].count; ++i) {
          if constexpr (collect_statistics) ++statistics->primitiveBoxTests;
          const int boxId           = primId[nodes[id].first + i];
          const Precision boxSafety = bvh.BoxSafety(boxId, localpoint);
          if (boxSafety >= safety) continue;
          if constexpr (collect_statistics) ++statistics->relevantBoxes;
          if (boxSafety > limit) {
            if constexpr (collect_statistics) {
              if (limit_shortcuts) ++*limit_shortcuts;
            }
            safety = vecgeom::Min(safety, boxSafety);
            continue;
          }
          if constexpr (collect_statistics) {
            if (candidate_calls) ++*candidate_calls;
          }
          const int prim       = bvh.GetDaughterId(boxId);
          const Precision dist = Navigator::CandidateSafetyToIn(bvh.GetRootId(), prim, localpoint);
          if (dist > min_valid_safety) safety = vecgeom::Min(safety, dist);
        }
      } else {
        const unsigned int childL = nodes[id].first;
        const unsigned int childR = childL + 1;
        const BVHReal_t safetyL   = bvh.NodeSafetySqr(childL, point);
        const BVHReal_t safetyR   = bvh.NodeSafetySqr(childR, point);
        const Precision limitSqr  = safety * safety;
        const bool traverseL      = static_cast<Precision>(safetyL) < limitSqr;
        const bool traverseR      = static_cast<Precision>(safetyR) < limitSqr;
        if (safetyR < safetyL) {
          if (traverseL) stack[stackSize++] = childL;
          if (traverseR) stack[stackSize++] = childR;
        } else {
          if (traverseR) stack[stackSize++] = childR;
          if (traverseL) stack[stackSize++] = childL;
        }
      }
    }
    return safety;
  }

protected:
  // Exact daughter operations used by the navigator traversal and specialized derived adapters.

  /*
   * @param[in] aLVIndex Global index of a LogicalVolume
   * @param[in] index Index within the list of daughters of the specified LogicalVolume
   * @returns The PlacedVolume defined by @p aLVIndex and @p index
   */
  VECCORE_ATT_HOST_DEVICE
  static VECGEOM_FORCE_INLINE Daughter GetPlacedVolume(int aLVIndex, int index)
  {
#ifdef VECCORE_CUDA_DEVICE_COMPILATION
    VECGEOM_VALIDATE(vecgeom::globaldevicegeomdata::gDeviceLogicalVolumes != nullptr,
                     << "Logical volumes not copied to device");
    return vecgeom::globaldevicegeomdata::gDeviceLogicalVolumes[aLVIndex].GetDaughters()[index];
#else
#ifndef VECCORE_CUDA
    return vecgeom::GeoManager::Instance().GetLogicalVolume(aLVIndex)->GetDaughters()[index];
#else
    // this is the case when we compile with nvcc for host side
    VECGEOM_VALIDATE(false, << "reached unimplement code");
    (void)index; // avoid unused parameter warning.
    (void)aLVIndex;
    return nullptr;
#endif
#endif
  }

  /*
   * @param[in] global_index Global index of a PlacedVolume
   * @returns The PlacedVolume with global index @p global_index
   */
  VECCORE_ATT_HOST_DEVICE
  static VECGEOM_FORCE_INLINE vecgeom::VPlacedVolume *GetPlacedVolume(int global_index)
  {
#ifdef VECCORE_CUDA_DEVICE_COMPILATION
    VECGEOM_VALIDATE(vecgeom::globaldevicegeomdata::gCompactPlacedVolBuffer != nullptr,
                     << "Placed volumes not copied to device");
    return &vecgeom::globaldevicegeomdata::gCompactPlacedVolBuffer[global_index];
#else
#ifndef VECCORE_CUDA
    return vecgeom::GeoManager::Instance().GetPlacedVolume(global_index);
#else
    // this is the case when we compile with nvcc for host side
    VECGEOM_VALIDATE(false, << "reached unimplement code");
    (void)global_index; // avoid unused parameter warning.
    return nullptr;
#endif
#endif
  }

  /*
   * @param[in] aLVIndex Global index of a LogicalVolume
   * @param[in] index Index within the list of daughters of the specified LogicalVolume
   * @param[in] localpoint Point in the local coordinates of the LV specified by @aLVIndex
   * @returns The safety to in to the PlacedVolume defined by @p aLVIndex and @p index for the point @p localpoint
   */
  VECCORE_ATT_HOST_DEVICE
  static Precision CandidateSafetyToIn(int aLVIndex, int index, Vector3D<Precision> localpoint)
  {
    return GetPlacedVolume(aLVIndex, index)->SafetyToIn(localpoint);
  };

  /*
   * @param[in] aLVIndex Global index of a LogicalVolume
   * @param[in] index Index within the list of daughters of the specified LogicalVolume
   * @param[in] localpoint Point in the local coordinates of the LV specified by @aLVIndex
   * @param[in] localdir Direction in the local coordinates of the LV specified by @aLVIndex
   * @param[in] step Maximum step length
   * @returns The distance to in to the PlacedVolume defined by @p aLVIndex and @p index for the point @p localpoint
   * and direction @p localdir
   */
  VECCORE_ATT_HOST_DEVICE
  static Precision CandidateDistanceToIn(int aLVIndex, int index, Vector3D<Precision> localpoint,
                                         Vector3D<Precision> localdir, Precision step)
  {
    Daughter vol = GetPlacedVolume(aLVIndex, index);
    return vol->DistanceToIn(localpoint, localdir, step);
  };

  /*
   * @param[in] aLVIndex Global index of a LogicalVolume
   * @param[in] index Index within the list of daughters of the specified LogicalVolume
   * @param[in] localpoint Point in the local coordinates of the LV specified by @aLVIndex, on or inside the daughter
   * @param[in] localdir Direction in the local coordinates of the LV specified by @aLVIndex
   * @returns Whether the surface normal of the PlacedVolume defined by @p aLVIndex and @p index at @p localpoint
   * has a positive projection on @p localdir, that is whether the ray is leaving it
   */
  VECCORE_ATT_HOST_DEVICE
  static bool CandidateLeaves(int aLVIndex, int index, Vector3D<Precision> localpoint, Vector3D<Precision> localdir)
  {
    Vector3D<Precision> normal;
    GetPlacedVolume(aLVIndex, index)->Normal(localpoint, normal);
    return normal.Dot(localdir) > 0.;
  }

  /*
   * @param[in] aLVIndex Global index of a LogicalVolume
   * @param[in] index Index within the list of daughters of the specified LogicalVolume
   * @param[in] localpoint Point in the local coordinates of the LV specified by @aLVIndex
   * @param[out] daughterlocalpoint Point in the local coordinates of the PlacedVolume defined by
   * @p aLVIndex and @p index
   * @returns Whether @localpoint falls within the PlacedVolume defined by @p aLVIndex and @p index
   */
  VECCORE_ATT_HOST_DEVICE
  static vecgeom::Inside_t CandidateInside(int aLVIndex, int index, Vector3D<Precision> const &localpoint,
                                           Vector3D<Precision> &daughterlocalpoint)
  {
    auto daughter      = GetPlacedVolume(aLVIndex, index);
    daughterlocalpoint = daughter->GetTransformation()->template Transform<Precision>(localpoint);
    return daughter->GetUnplacedVolume()->Inside(daughterlocalpoint);
  };

  /*
   * @param[in] aLVIndex Global index of a LogicalVolume
   * @param[in] index Index within the list of daughters of the specified LogicalVolume
   * @param[in] localpoint Point in the local coordinates of the LV specified by @aLVIndex
   * @param[out] daughterlocalpoint Point in the local coordinates of the PlacedVolume defined by
   * @p aLVIndex and @p index
   * @returns Whether @localpoint falls within the PlacedVolume defined by @p aLVIndex and @p index
   */
  VECCORE_ATT_HOST_DEVICE
  static bool CandidateContains(int aLVIndex, int index, Vector3D<Precision> const &localpoint,
                                Vector3D<Precision> &daughterlocalpoint)
  {
    auto inside = CandidateInside(aLVIndex, index, localpoint, daughterlocalpoint);
    return inside != EnumInside::kOutside;
  };

  /*
   * @param[in] aLVIndex Global index of a LogicalVolume
   * @param[in] index Index within the list of daughters of the specified LogicalVolume
   * @param[in] localpoint Point in the local coordinates of the LV specified by @aLVIndex
   * @param[in] localdir Direction in the local coordinates of the LV specified by @aLVIndex
   * @returns The distance to in to the Bounding Box of the PlacedVolume defined by @p aLVIndex
   * and @p index for the point @p localpoint and direction @p localdir
   */
  VECCORE_ATT_HOST_DEVICE
  static Precision CandidateApproachSolid(int aLVIndex, int index, Vector3D<Precision> localpoint,
                                          Vector3D<Precision> localdir)
  {
    auto vol                            = GetPlacedVolume(aLVIndex, index);
    vecgeom::Transformation3D const *tr = vol->GetTransformation();
    Vector3D<Precision> pv_localpoint   = tr->Transform(localpoint);
    Vector3D<Precision> pv_localdir     = tr->TransformDirection(localdir);
    Vector3D<Precision> pv_invlocaldir(1.0 / vecgeom::NonZero(pv_localdir[0]), 1.0 / vecgeom::NonZero(pv_localdir[1]),
                                       1.0 / vecgeom::NonZero(pv_localdir[2]));
    return vol->GetUnplacedVolume()->ApproachSolid(pv_localpoint, pv_invlocaldir);
  };

  /*
   * Used by the BVH to determine if it needs to skip checking a placed volume. The global index of the volume
   * defined by @p aLVIndex and @p index can only be accessed from the navigator
   * @param[in] aLVIndex Global index of a LogicalVolume
   * @param[in] index Index within the list of daughters of the specified LogicalVolume
   * @param[in] global_id Global id of a PLacedVolume
   * @returns Whether the global id of the PlacedVolume defined by @p aLVIndex and @p index is the same as @p global_id
   */
  VECCORE_ATT_HOST_DEVICE
  static VECGEOM_FORCE_INLINE bool SkipItem(int aLVIndex, int index, long const global_id)
  {
    return (global_id == GetPlacedVolume(aLVIndex, index)->id());
  }

public:
  /** @name Internal component interface
   *  These traversal operations are shared with VecGeom safety and surface-navigation components.
   *  @{ */

  /**
   * @brief Find the nearest daughter entry along a bounded IndexedBVH ray query.
   *
   * Primitive-box ids are mapped before the exact query. The exact query starts
   * at the box entry distance; mapped boxes clamp a negative residual distance
   * because a split box may begin inside its daughter.
   *
   * @tparam Navigator Static adapter providing exact daughter operations.
   * @tparam BVHReal_t Scalar type stored by the hierarchy.
   * @tparam Real_i Scalar query precision.
   * @param bvh Box hierarchy for the mother volume.
   * @param localpoint Ray origin in mother-local coordinates.
   * @param localdir Ray direction in mother-local coordinates.
   * @param[in,out] step Maximum step on input; nearest accepted entry on output.
   * @param last_exited_id Placed-volume id excluded from immediate re-entry, or
   *                       a negative value when no exclusion applies.
   * @param[out] hitcandidate_index Daughter index of the nearest accepted hit;
   *                                unchanged when no closer daughter is found.
   */
  template <typename Navigator, typename BVHReal_t, typename Real_i>
  VECCORE_ATT_HOST_DEVICE static void CheckDaughterIntersections(const vecgeom::IndexedBVH<BVHReal_t> &bvh,
                                                                 const Vector3D<Real_i> &localpoint,
                                                                 const Vector3D<Real_i> &localdir, Real_i &step,
                                                                 long const last_exited_id, long &hitcandidate_index)
  {
    CheckDaughterIntersectionsImpl<false, Navigator>(bvh, localpoint, localdir, step, last_exited_id,
                                                     hitcandidate_index);
  }

  /**
   * @brief CheckDaughterIntersections for Geant4: at a non-positive distance, refuses @p blocked_id and any
   * daughter the ray is leaving, as SimpleABBoxNavigator does.
   * @param blocked_id Placed-volume id refused at a non-positive distance, or -1 for none.
   */
  template <typename Navigator, typename BVHReal_t, typename Real_i>
  VECCORE_ATT_HOST_DEVICE static void CheckDaughterIntersectionsBlocked(const vecgeom::IndexedBVH<BVHReal_t> &bvh,
                                                                        const Vector3D<Real_i> &localpoint,
                                                                        const Vector3D<Real_i> &localdir, Real_i &step,
                                                                        long const blocked_id, long &hitcandidate_index)
  {
    CheckDaughterIntersectionsImpl<false, Navigator, true>(bvh, localpoint, localdir, step, blocked_id,
                                                           hitcandidate_index);
  }

  /**
   * @brief Run `CheckDaughterIntersections()` while collecting traversal work.
   *
   * This diagnostic query preserves the candidate semantics of the normal
   * query but is not used by navigation.
   *
   * @tparam Navigator Static adapter providing exact daughter operations.
   * @tparam BVHReal_t Scalar type stored by the hierarchy.
   * @tparam Real_i Scalar query precision.
   * @param bvh Box hierarchy for the mother volume.
   * @param localpoint Ray origin in mother-local coordinates.
   * @param localdir Ray direction in mother-local coordinates.
   * @param[in,out] step Maximum step on input; nearest accepted entry on output.
   * @param last_exited_id Placed-volume id excluded from immediate re-entry.
   * @param[out] hitcandidate_index Daughter index of the nearest accepted hit.
   * @param[out] total_visited_children Primitive AABB tests.
   * @param[out] total_visited_leaves Visited leaf nodes.
   * @param[out] total_cut_nodes Root-node rejections.
   * @param[out] total_scheduled_nodes Intersected child nodes.
   * @param[out] total_intersected_children Intersected primitive AABBs, when requested.
   * @param[out] total_candidate_calls Exact daughter calls, when requested.
   * @param[out] total_deferred_nodes Deferred far children, when requested.
   * @param[out] total_prunable_deferred_nodes Deferred children beyond a reduced step, when requested.
   * @param[out] maximum_stack_depth Maximum deferred-node count, when requested.
   */
  template <typename Navigator, typename BVHReal_t, typename Real_i>
  VECCORE_ATT_HOST_DEVICE static void CheckDaughterIntersectionsBenchmark(
      const vecgeom::IndexedBVH<BVHReal_t> &bvh, const Vector3D<Real_i> &localpoint, const Vector3D<Real_i> &localdir,
      Real_i &step, long const last_exited_id, long &hitcandidate_index, long &total_visited_children,
      long &total_visited_leaves, long &total_cut_nodes, long &total_scheduled_nodes,
      long *total_intersected_children = nullptr, long *total_candidate_calls = nullptr,
      long *total_deferred_nodes = nullptr, long *total_prunable_deferred_nodes = nullptr,
      long *maximum_stack_depth = nullptr)
  {
    IntersectStatistics statistics;
    if (total_candidate_calls) *total_candidate_calls = 0;
    CheckDaughterIntersectionsImpl<true, Navigator>(bvh, localpoint, localdir, step, last_exited_id, hitcandidate_index,
                                                    &statistics, total_candidate_calls);

    total_visited_children = statistics.primitiveBoxTests;
    total_visited_leaves   = statistics.visitedLeaves;
    total_cut_nodes        = statistics.cutNodes;
    total_scheduled_nodes  = statistics.scheduledNodes;
    if (total_intersected_children) *total_intersected_children = statistics.intersectedBoxes;
    if (total_deferred_nodes) *total_deferred_nodes = statistics.deferredNodes;
    if (total_prunable_deferred_nodes) *total_prunable_deferred_nodes = statistics.prunableDeferredNodes;
    if (maximum_stack_depth) *maximum_stack_depth = statistics.maximumStackDepth;
  }

  /**
   * @brief Reduce @p step to the nearest approached daughter solid.
   *
   * This method traverses primitive boxes near-first, maps each candidate, and
   * applies the reduced step immediately to prune the remaining traversal.
   *
   * @tparam Navigator Static adapter providing exact daughter operations.
   * @tparam BVHReal_t Scalar type stored by the hierarchy.
   * @param bvh Box hierarchy for the mother volume.
   * @param localpoint Ray origin in mother-local coordinates.
   * @param localdir Ray direction in mother-local coordinates.
   * @param[in,out] step Maximum approach on input; nearest accepted approach on output.
   * @param last_exited_id Placed-volume id excluded at non-positive distance.
   */
  template <typename Navigator, typename BVHReal_t>
  VECCORE_ATT_HOST_DEVICE static void ApproachNextDaughter(const vecgeom::IndexedBVH<BVHReal_t> &bvh,
                                                           Vector3D<Precision> const &localpoint,
                                                           Vector3D<Precision> const &localdir, Precision &step,
                                                           long const last_exited_id)
  {
    unsigned int stack[vecgeom::IndexedBVH<BVHReal_t>::kTraversalStackSize];
    int stackSize = 0;
    const Vector3D<BVHReal_t> point(static_cast<BVHReal_t>(localpoint[0]), static_cast<BVHReal_t>(localpoint[1]),
                                    static_cast<BVHReal_t>(localpoint[2]));
    const Vector3D<BVHReal_t> invDirection(static_cast<BVHReal_t>(1.) / vecgeom::NonZero(localdir[0]),
                                           static_cast<BVHReal_t>(1.) / vecgeom::NonZero(localdir[1]),
                                           static_cast<BVHReal_t>(1.) / vecgeom::NonZero(localdir[2]));
    const unsigned int octant = (invDirection[0] < BVHReal_t(0.) ? 1u : 0u) |
                                (invDirection[1] < BVHReal_t(0.) ? 2u : 0u) |
                                (invDirection[2] < BVHReal_t(0.) ? 4u : 0u);
    BVHReal_t bvhStep         = static_cast<BVHReal_t>(step);
    BVHReal_t approach;
    if (!bvh.IntersectNode(0, point, invDirection, octant, bvhStep, approach)) return;

    auto const *nodes  = bvh.GetNodes();
    auto const *primId = bvh.GetPrimId();
    unsigned int top   = 0;
    for (;;) {
      while (!nodes[top].IsLeaf()) {
        const unsigned int childL = nodes[top].first;
        const unsigned int childR = childL + 1;
        BVHReal_t approachL, approachR;
        const bool hitL = bvh.IntersectNode(childL, point, invDirection, octant, bvhStep, approachL);
        const bool hitR = bvh.IntersectNode(childR, point, invDirection, octant, bvhStep, approachR);
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
        if (!bvh.IntersectBox(boxId, point, invDirection, bvhStep, approach)) continue;
        const int prim  = bvh.GetDaughterId(boxId);
        const auto dist = Navigator::CandidateApproachSolid(bvh.GetRootId(), prim, localpoint, localdir);
        if (dist < step && !(dist <= Precision(0.) && Navigator::SkipItem(bvh.GetRootId(), prim, last_exited_id))) {
          step    = dist;
          bvhStep = static_cast<BVHReal_t>(dist);
        }
      }

    pop:
      if (stackSize == 0) return;
      top = stack[--stackSize];
    }
  }

  /**
   * @brief Reduce a Safety bound against mapped primitive candidates.
   *
   * A primitive whose conservative AABB Safety exceeds @p limit can reduce the
   * result without an exact item query. Other accepted boxes are mapped before
   * the exact query.
   *
   * @tparam Navigator Static adapter providing exact item Safety operations.
   * @tparam BVHReal_t Scalar type stored by the hierarchy.
   * @param bvh Box hierarchy for the owner.
   * @param localpoint Point in owner-local coordinates.
   * @param safety Initial upper bound.
   * @param limit Accuracy limit for exact item calls.
   * @param min_valid_safety Ignore exact values at or below this value.
   * @return Minimum accepted item Safety or the input @p safety.
   */
  template <typename Navigator, typename BVHReal_t>
  VECCORE_ATT_HOST_DEVICE static Precision ComputeBVHSafety(const vecgeom::IndexedBVH<BVHReal_t> &bvh,
                                                            Vector3D<Precision> const &localpoint, Precision safety,
                                                            Precision limit            = InfinityLength<Precision>(),
                                                            Precision min_valid_safety = -InfinityLength<Precision>())
  {
    return ComputeBVHSafetyImpl<false, Navigator>(bvh, localpoint, safety, limit, min_valid_safety);
  }

  /**
   * @brief Run `ComputeBVHSafety()` while collecting traversal work.
   *
   * This diagnostic query uses the tree path on the CPU so its counters describe
   * actual nodes and leaves; it is not used by navigation.
   *
   * @tparam Navigator Static adapter providing exact item Safety operations.
   * @tparam BVHReal_t Scalar type stored by the hierarchy.
   * @param bvh Box hierarchy for the owner.
   * @param localpoint Point in owner-local coordinates.
   * @param safety Initial upper bound.
   * @param limit Accuracy limit for exact item calls.
   * @param min_valid_safety Ignore exact values at or below this value.
   * @param[out] visited_nodes Visited hierarchy nodes.
   * @param[out] visited_leaves Visited leaf nodes.
   * @param[out] primitive_box_tests Primitive AABB Safety evaluations.
   * @param[out] relevant_boxes Primitive AABBs closer than the current bound.
   * @param[out] limit_shortcuts AABB results accepted above @p limit.
   * @param[out] candidate_calls Exact item Safety calls.
   * @param[out] cut_nodes Nodes rejected by the current bound.
   * @return Reduced Safety.
   */
  template <typename Navigator, typename BVHReal_t>
  VECCORE_ATT_HOST_DEVICE static Precision ComputeSafetyBenchmark(
      const vecgeom::IndexedBVH<BVHReal_t> &bvh, Vector3D<Precision> const &localpoint, Precision safety,
      Precision limit, Precision min_valid_safety, long &visited_nodes, long &visited_leaves, long &primitive_box_tests,
      long &relevant_boxes, long &limit_shortcuts, long &candidate_calls, long &cut_nodes)
  {
    SafetyStatistics statistics;
    limit_shortcuts = 0;
    candidate_calls = 0;
    safety        = ComputeBVHSafetyImpl<true, Navigator>(bvh, localpoint, safety, limit, min_valid_safety, &statistics,
                                                          &limit_shortcuts, &candidate_calls);
    visited_nodes = statistics.visitedNodes;
    visited_leaves      = statistics.visitedLeaves;
    primitive_box_tests = statistics.primitiveBoxTests;
    relevant_boxes      = statistics.relevantBoxes;
    cut_nodes           = statistics.cutNodes;
    return safety;
  }

  /**
   * @brief Find the first non-excluded daughter containing @p localpoint.
   *
   * @tparam Navigator Static adapter providing exact daughter operations and ids.
   * @tparam BVHReal_t Scalar type stored by the hierarchy.
   * @param bvh Box hierarchy for the mother volume.
   * @param exclude_item_id Placed-volume id to skip, or a negative value.
   * @param localpoint Point in mother-local coordinates.
   * @param[out] container_id Placed-volume id of the containing daughter.
   * @param[out] daughterlocalpoint Point transformed into that daughter.
   * @return `kInside` or `kSurface` for the first accepted daughter, otherwise `kOutside`.
   */
  template <typename Navigator, typename BVHReal_t>
  VECCORE_ATT_HOST_DEVICE static vecgeom::Inside_t LevelInside(const vecgeom::IndexedBVH<BVHReal_t> &bvh,
                                                               long const exclude_item_id,
                                                               Vector3D<Precision> const &localpoint,
                                                               long &container_id,
                                                               Vector3D<Precision> &daughterlocalpoint)
  {
    vecgeom::Inside_t result = kOutside;
    const Vector3D<BVHReal_t> point(static_cast<BVHReal_t>(localpoint[0]), static_cast<BVHReal_t>(localpoint[1]),
                                    static_cast<BVHReal_t>(localpoint[2]));
    if (!bvh.NodeContains(0, point)) return result;

    unsigned int stack[vecgeom::IndexedBVH<BVHReal_t>::kTraversalStackSize];
    int stackSize      = 0;
    auto const *nodes  = bvh.GetNodes();
    auto const *primId = bvh.GetPrimId();
    unsigned int top   = 0;
    for (;;) {
      while (!nodes[top].IsLeaf()) {
        const unsigned int childL = nodes[top].first;
        const unsigned int childR = childL + 1;
        const bool inL            = bvh.NodeContains(childL, point);
        const bool inR            = bvh.NodeContains(childR, point);
        if (inL) {
          if (inR) stack[stackSize++] = childR;
          top = childL;
        } else if (inR) {
          top = childR;
        } else {
          goto pop;
        }
      }

      for (int i = 0; i < nodes[top].count; ++i) {
        const int boxId = primId[nodes[top].first + i];
        if (!bvh.BoxContains(boxId, point)) continue;
        const int prim = bvh.GetDaughterId(boxId);
        if (Navigator::SkipItem(bvh.GetRootId(), prim, exclude_item_id)) continue;
        const auto inside = Navigator::CandidateInside(bvh.GetRootId(), prim, localpoint, daughterlocalpoint);
        if (inside == kOutside) continue;
        container_id = Navigator::ItemId(bvh.GetRootId(), prim);
        return inside;
      }

    pop:
      if (stackSize == 0) return result;
      top = stack[--stackSize];
    }
  }

  /** @} */

protected:
  /** @name Navigation diagnostics
   *  Instrumented entry points used by validation and benchmarking.
   *  @{ */

  /** Test adapter exposing one IndexedBVH daughter-intersection query. */
  VECCORE_ATT_HOST_DEVICE
  static long TestBVHCheckDaughterIntersections(const vecgeom::IndexedBVH<Real_t> &bvh, Vector3D<Precision> &localpoint,
                                                Vector3D<Precision> &localdir, Precision &bvhstep)
  {
    long hitcandidate_index = -1;
    long last_exited_id     = -1;
    CheckDaughterIntersections<Self_t>(bvh, localpoint, localdir, bvhstep, last_exited_id, hitcandidate_index);
    return hitcandidate_index;
  }

  VECCORE_ATT_HOST_DEVICE
  static long TestCheckDaughterIntersections(int logicalVolumeId, Vector3D<Precision> const &localpoint,
                                             Vector3D<Precision> const &localdir, Precision &bvhstep)
  {
    long hitcandidate_index = -1;
    auto const *bvh         = GetBVH(logicalVolumeId);
    CheckDaughterIntersections<Self_t>(*bvh, localpoint, localdir, bvhstep, -1, hitcandidate_index);
    return hitcandidate_index;
  }

  /** @} */

  /*
   * @param[in] aLVIndex Global index of a LogicalVolume
   * @param[in] index Index within the list of daughters of the specified LogicalVolume
   * @returns The global id of the PlacedVolume defined by @p aLVIndex and @p index
   */
  VECCORE_ATT_HOST_DEVICE
  static uint ItemId(int aLVIndex, int index) { return GetPlacedVolume(aLVIndex, index)->id(); }

private:
  /** @name Point-location implementation
   *  These bodies resolve point-location queries using the internal BVH traversal.
   *  @{ */

  /**
   * @brief Locate the deepest daughter volume that contains a point, starting from a given placed volume.
   *
   * ## Coordinate frames
   * - **Input point** @p point is in the **parent frame of @p vol** (i.e. the same frame expected by
   *   `VPlacedVolume::Inside` for @p vol).
   * - The NavigationState @p path should point to the **parent frame of @p vol**.
   * @param[in]  vol
   *   The starting placed volume (the mother at which the descent begins). Must be non-null if @p top is true.
   * @param[in]  point
   *   Query point **in the parent frame of @p vol**.
   * @param[in,out] path
   *   Navigation state to be filled. Must start at the  **parent frame of @p vol**.
   * @param[in]  top
   *   If true, the function first validates that the point is inside @p vol. If false, a nullptr is returned
   * @param[in]  exclude
   *   Optional placed volume to exclude once during the descent (useful to avoid immediately
   *   re-entering a volume you just exited). The exclusion is applied only to the first BVH
   *   query and then cleared.
   * @return
   *   The deepest placed volume that contains the point, or `nullptr` if @p top is true and
   *   `vol->Inside(point)` reports outside. When the descent stops at @p vol (i.e. no daughter
   *   contains the point), returns @p vol.
   */
  VECCORE_ATT_HOST_DEVICE
  static Daughter LocatePointInImpl(vecgeom::VPlacedVolume const *vol, Vector3D<Precision> const &point,
                                    vecgeom::NavigationState &path, bool top,
                                    vecgeom::VPlacedVolume const *exclude = nullptr)
  {

    // optional check whether the point is inside the provided placed volume vol
    if (top) {
      // Must check the provided volume
      VECGEOM_ASSERT(vol != nullptr);
      auto inside = vol->Inside(point);
      if (inside == kOutside) return nullptr;
      // Set the boundary state to the path
      if (inside == kSurface) path.SetBoundaryState(true);
    }

    path.Push(vol);

    // transform the point into the reference frame of the `vol`
    Vector3D<Precision> currentpoint(vol->GetTransformation()->Transform<Precision>(point));
    Vector3D<Precision> daughterlocalpoint;
    long exclude_id = -1;
    long vol_id     = -1;

    for (auto v = vol; v->GetDaughters().size() > 0;) {
      auto bvh = GetBVH(v->GetLogicalVolume()->id());

      exclude_id = -1;
      if (exclude != nullptr) {
        exclude_id = exclude->id();
      }
      vol_id = -1;

      const auto inside = LevelInside<Self_t>(*bvh, exclude_id, currentpoint, vol_id, daughterlocalpoint);
      if (inside == kOutside) break;
      if (inside == kSurface) path.SetBoundaryState(true);

      currentpoint = daughterlocalpoint;
      // Update the current volume v
      v = GetPlacedVolume(vol_id);
      path.Push(v);

      // Only exclude the placed volume once since we could enter it again via a
      // different volume history.
      exclude = nullptr;
    }

    return path.Top();
  }

  /**
   * @brief Find new final NavigationState after a NavigationState has just been left (i.e., after exiting its corresponding top volume)
   * The function will climb up the NavigationStates until the point is inside, and then descend again to find the
   * deepest volume containing the point. The just exited volume will be skipped
   * @param[in]  localpoint   Point in the **local frame of the navigation state (path.Top())**
   * @param[in,out] path      Navigation state of the volume that was just left; updated to the new deepest containing path.
   * @return The deepest placed volume that contains the point (after re-location), or nullptr if no mother exists.
   */
  VECCORE_ATT_HOST_DEVICE
  static Daughter RelocatePointImpl(Vector3D<Precision> const &localpoint, vecgeom::NavigationState &path)
  {
    vecgeom::VPlacedVolume const *currentmother = path.Top();
    Daughter skip                               = nullptr;
    Vector3D<Precision> transformed             = localpoint;
    // continue to climb up the Navigation tree until a volume is found which contains the volume
    do {
      skip = currentmother;
      path.Pop();
      transformed   = currentmother->GetTransformation()->InverseTransform(transformed);
      currentmother = path.Top();
    } while (currentmother && (currentmother->IsAssembly() || !currentmother->UnplacedContains(transformed)));

    // first volume up the hierarchy is found that contains the point. Now descend to find the deepest state that
    // contains the point
    if (currentmother) {
      return LocatePointInNavState(transformed, path, false, skip);
    }
    return currentmother;
  }

  /**
   * @brief Update NavState path to the deepest state that contains given point starting from the given path.
   * @param[in]  localpoint
   *   localpoint **in the reference frame of @p path.**.
   * @param[in,out] path
   *   Starting NavigationState; output of the final NavigationState
   * @param[in]  top
   *   If true, the function first validates that the point is inside @p path. If false, a nullptr is returned
   * @param[in]  exclude
   *   Optional placed volume to exclude once during the descent
   * @return
   *   The final top-volume of the NavState after finding the deepest NavState that contains the point
   */
  VECCORE_ATT_HOST_DEVICE
  static Daughter LocatePointInNavStateImpl(Vector3D<Precision> const &localpoint, vecgeom::NavigationState &path,
                                            bool top, vecgeom::VPlacedVolume const *exclude = nullptr)
  {

    VECGEOM_ASSERT(path.Top() != nullptr);

    auto currentLogical = path.Top()->GetLogicalVolume();
    VECGEOM_ASSERT(currentLogical != nullptr);

    // optional check if point is inside current NavState's top volume
    if (top) {
      auto inside = currentLogical->GetUnplacedVolume()->Inside(localpoint);
      if (inside == kOutside) return nullptr;
      // Set the boundary state to the path
      if (inside == kSurface) path.SetBoundaryState(true);
    }

    Vector3D<Precision> currentpoint(localpoint);
    Vector3D<Precision> daughterlocalpoint;
    // exclude->id() returns unsigned int, but LevelInside takes a long as an input
    long exclude_id = exclude ? static_cast<long>(exclude->id()) : -1;
    long pvol_id    = -1;

    while (currentLogical->GetDaughters().size() > 0) {
      auto *bvh = GetBVH(currentLogical->id());

      const auto inside = LevelInside<Self_t>(*bvh, exclude_id, currentpoint, pvol_id, daughterlocalpoint);

      if (inside == kOutside) break; // no containing daughter here
      if (inside == kSurface) path.SetBoundaryState(true);

      auto *placed = GetPlacedVolume(pvol_id);
      path.Push(placed);

      // Prepare for next level: switch to daughter's logical volume + daughter-local point
      currentpoint   = daughterlocalpoint;
      currentLogical = placed->GetLogicalVolume();

      // Only exclude the placed volume once since we could enter it again via a
      // different volume history.
      exclude_id = -1;
      exclude    = nullptr;
    }

    return path.Top();
  }

  /** @} */

  // Computes a step in the current volume from the localpoint into localdir,
  // taking step_limit into account. If a volume is hit, the function calls
  // out_state.SetBoundaryState(true) and hitcandidate is set to the hit
  // daughter volume, or kept unchanged if the current volume is left.
  VECCORE_ATT_HOST_DEVICE
  static Precision ComputeStepAndHit(Vector3D<Precision> const &localpoint, Vector3D<Precision> const &localdir,
                                     Precision step_limit, vecgeom::NavigationState const &in_state,
                                     vecgeom::NavigationState &out_state, Daughter &hitcandidate)
  {
    in_state.CopyTo(&out_state);
    // Just return unchanged current state if the step limit is null or invalid
    if (step_limit <= 0) return 0.;

    Precision step = step_limit;
    Daughter pvol  = in_state.Top();
    // Daughter last_exited = in_state.GetLastExited();
    long hitcandidate_index = -1;
    long last_exited_id     = -1;

    // need to calc DistanceToOut first
    step = pvol->DistanceToOut(localpoint, localdir, step_limit);
    // This should never happen as DistanceToOut should never return infinity (TBC)
    if (step == vecgeom::kInfLength) step = 0.;
    // Most solids ignore step_limit
    step = Min(step, step_limit);
    // Boundary very close or already outside
    if (step < kTolerance) {
      step = Max(Precision(0.), step);
      // Even if the point is outside the volume, we have to set the boundary flag to force relocation
      out_state.SetBoundaryState(true);
      return step;
    }

    // Now distance to children
    if (pvol->GetDaughters().size() > 0) {
      auto bvh = GetBVH(pvol->GetLogicalVolume()->id());

      hitcandidate_index = -1;
      // id is an uint, however we use a long in order to be able to fit the full uint range, and -1 in case there is no
      // last exited volume in the navigation state.
      last_exited_id = -1;
      // if (last_exited != nullptr) last_exited_id = last_exited->id();

      CheckDaughterIntersections<Self_t>(*bvh, localpoint, localdir, step, last_exited_id, hitcandidate_index);

      if (hitcandidate_index >= 0) {
        // A child was hit within the step_limit
        step         = Max(Precision(0.), step);
        hitcandidate = pvol->GetLogicalVolume()->GetDaughters()[hitcandidate_index];
        out_state.SetBoundaryState(true);
        return step;
      }
    }

    // Is geometry further away than physics step?
    if (step >= step_limit) {
      // Then this is a phyics step and we don't need to do anything.
      out_state.SetBoundaryState(false);
      return step_limit;
    }

    // Otherwise it is a geometry step and we push the point to the boundary.
    out_state.SetBoundaryState(true);
    return step;
  }

  // Computes a step in the current volume from the localpoint into localdir,
  // until the next daughter bounding box, taking step_limit into account.
  VECCORE_ATT_HOST_DEVICE
  static Precision ApproachNextVolume(Vector3D<Precision> const &localpoint, Vector3D<Precision> const &localdir,
                                      Precision step_limit, vecgeom::NavigationState const &in_state)
  {
    Precision step = step_limit;
    Daughter pvol  = in_state.Top();
    // Daughter last_exited = in_state.GetLastExited();

    if (pvol->GetDaughters().size() > 0) {
      auto bvh = GetBVH(pvol->GetLogicalVolume()->id());

      // id is an uint, however we use a long in order to be able to fit the full uint range, and -1 in case there is no
      // last exited volume in the navigation state.
      long last_exited_id = -1;
      // if (last_exited != nullptr) last_exited_id = last_exited->id();

      ApproachNextDaughter<Self_t>(*bvh, localpoint, localdir, step, last_exited_id);
      // Make sure we don't "step" on next boundary
      step -= 10 * vecgeom::kTolerance;
    }

    if (step == vecgeom::kInfLength && step_limit > 0) return 0;

    // Is geometry further away than physics step?
    if (step > step_limit) {
      // Then this is a phyics step and we don't need to do anything.
      return step_limit;
    }

    if (step < 0) {
      step = 0;
    }

    return step;
  }

  /** @name Tracking implementation
   *  These bodies transform global queries and compose the internal traversal helpers.
   *  @{ */

  /**
   * @brief Compute isotropic Safety from a global point in the current state.
   * @param globalpoint Point in global coordinates.
   * @param state Current navigation state.
   * @param limit Accuracy limit; values above it may be conservative.
   * @return Safety to the nearest boundary, or infinity when @p state is outside.
   */
  VECCORE_ATT_HOST_DEVICE
  static Precision ComputeSafetyImpl(Vector3D<Precision> const &globalpoint, vecgeom::NavigationState const &state,
                                     Precision limit)
  {
    Daughter pvol = state.Top();
    if (pvol == nullptr) return kInfLength;
    vecgeom::Transformation3D m;
    state.TopMatrix(m);
    Vector3D<Precision> localpoint = m.Transform(globalpoint);

    // need to calc DistanceToOut first
    Precision safety = pvol->SafetyToOut(localpoint);
    limit            = Min(safety, limit);

    if (safety > 0 && pvol->GetDaughters().size() > 0) {
      auto bvh = GetBVH(pvol->GetLogicalVolume()->id());
      safety   = ComputeBVHSafety<Self_t>(*bvh, localpoint, safety, limit);
    }

    return safety;
  }

  /**
   * @brief Compute a bounded step and fully relocate the output state.
   * @param globalpoint Point in global coordinates and inside @p in_state.
   * @param globaldir Unit direction in global coordinates.
   * @param step_limit Maximum requested step.
   * @param in_state Current navigation state.
   * @param[out] out_state State at the endpoint, including transitive boundary relocation.
   * @return Geometry-limited or requested step length.
   */
  VECCORE_ATT_HOST_DEVICE
  static Precision ComputeStepAndPropagatedStateImpl(Vector3D<Precision> const &globalpoint,
                                                     Vector3D<Precision> const &globaldir, Precision step_limit,
                                                     vecgeom::NavigationState const &in_state,
                                                     vecgeom::NavigationState &out_state)
  {
    if (in_state.Top() == nullptr) return kInfLength;
    const Precision push = in_state.IsOnBoundary() ? kBoundaryPush : 0.;

    if (step_limit < push) {
      // Go as far as the step limit says, assuming there is no boundary.
      // TODO: Does this make sense?
      in_state.CopyTo(&out_state);
      out_state.SetBoundaryState(false);
      return step_limit;
    }
    step_limit -= push;

    // calculate local point/dir from global point/dir
    Vector3D<Precision> localpoint;
    Vector3D<Precision> localdir;
    // Impl::DoGlobalToLocalTransformation(in_state, globalpoint, globaldir, localpoint, localdir);
    vecgeom::Transformation3D m;
    in_state.TopMatrix(m);
    localpoint = m.Transform(globalpoint);
    localdir   = m.TransformDirection(globaldir);

    Daughter hitcandidate = nullptr;
    // Avoid computing the distance from boundary by pushing the point
    Precision step =
        ComputeStepAndHit(localpoint + push * localdir, localdir, step_limit, in_state, out_state, hitcandidate);
    step += push;

    if (out_state.IsOnBoundary()) {
      // Relocate the point after the step to refine out_state.
      localpoint += (step + kBoundaryPush) * localdir;

      if (!hitcandidate) {
        // We didn't hit a daughter but instead we're exiting the current volume.
        RelocatePoint(localpoint, out_state);
      } else {
        // Otherwise check if we're directly entering other daughters transitively.
        LocatePointIn(hitcandidate, localpoint, out_state, false);
      }

      if (out_state.Top() != nullptr) {
        while (out_state.Top()->IsAssembly() || out_state.HasSamePathAsOther(in_state)) {
          out_state.Pop();
        }
        VECGEOM_ASSERT(!out_state.Top()->GetLogicalVolume()->GetUnplacedVolume()->IsAssembly());
      }
    }

    return step;
  }

  /**
   * @brief Compute a bounded step and update the state to the first crossed boundary.
   *
   * The output state enters the hit daughter or leaves the current volume, but
   * does not resolve additional coincident boundaries. Call
   * `RelocateToNextVolume()` after moving the track to the returned endpoint.
   *
   * @param globalpoint Point in global coordinates and inside @p in_state.
   * @param globaldir Unit direction in global coordinates.
   * @param step_limit Maximum requested step.
   * @param in_state Current navigation state.
   * @param[out] out_state State at the first crossed boundary, or a copy of @p in_state for a physics step.
   * @return Geometry-limited or requested step length.
   */
  VECCORE_ATT_HOST_DEVICE
  static Precision ComputeStepAndNextVolumeImpl(Vector3D<Precision> const &globalpoint,
                                                Vector3D<Precision> const &globaldir, Precision step_limit,
                                                vecgeom::NavigationState const &in_state,
                                                vecgeom::NavigationState &out_state)
  {
    if (in_state.Top() == nullptr) return kInfLength;
    const Precision push = in_state.IsOnBoundary() ? kBoundaryPush : 0.;

    if (step_limit < push) {
      // Go as far as the step limit says, assuming there is no boundary.
      // TODO: Does this make sense?
      in_state.CopyTo(&out_state);
      if (step_limit > kTolerance) out_state.SetBoundaryState(false);
      return step_limit;
    }
    step_limit -= push;

    // calculate local point/dir from global point/dir
    Vector3D<Precision> localpoint;
    Vector3D<Precision> localdir;
    // Impl::DoGlobalToLocalTransformation(in_state, globalpoint, globaldir, localpoint, localdir);
    vecgeom::Transformation3D m;
    in_state.TopMatrix(m);
    localpoint = m.Transform(globalpoint);
    localdir   = m.TransformDirection(globaldir);

    Daughter hitcandidate = nullptr;
    // Avoid computing the distance from boundary by pushing the point
    Precision step =
        ComputeStepAndHit(localpoint + push * localdir, localdir, step_limit, in_state, out_state, hitcandidate);
    step += (step > 0.) * push;

    if (out_state.IsOnBoundary()) {
      if (!hitcandidate) {
        vecgeom::VPlacedVolume const *currentmother = out_state.Top();
        Vector3D<Precision> transformed             = localpoint;
        // Push the point inside the next volume.
        transformed += (step + kBoundaryPush) * localdir;
        do {
          out_state.SetLastExited();
          out_state.Pop();
          transformed   = currentmother->GetTransformation()->InverseTransform(transformed);
          currentmother = out_state.Top();
        } while (currentmother &&
                 (currentmother->IsAssembly() || currentmother->GetUnplacedVolume()->Inside(transformed) != kInside));
      } else {
        out_state.Push(hitcandidate);
      }
    }

    return step;
  }

  /**
   * @brief Compute a bounded approach step without changing the navigation state.
   * @param globalpoint Point in global coordinates and inside @p in_state.
   * @param globaldir Unit direction in global coordinates.
   * @param step_limit Maximum requested step.
   * @param in_state Current navigation state.
   * @return Step limited by the next daughter approach or by @p step_limit.
   */
  VECCORE_ATT_HOST_DEVICE
  static Precision ComputeStepToApproachNextVolumeImpl(Vector3D<Precision> const &globalpoint,
                                                       Vector3D<Precision> const &globaldir, Precision step_limit,
                                                       vecgeom::NavigationState const &in_state)
  {
    if (in_state.Top() == nullptr) return kInfLength;
    // calculate local point/dir from global point/dir
    Vector3D<Precision> localpoint;
    Vector3D<Precision> localdir;
    // Impl::DoGlobalToLocalTransformation(in_state, globalpoint, globaldir, localpoint, localdir);
    vecgeom::Transformation3D m;
    in_state.TopMatrix(m);
    localpoint = m.Transform(globalpoint);
    localdir   = m.TransformDirection(globaldir);

    Precision step = ApproachNextVolume(localpoint, localdir, step_limit, in_state);

    return step;
  }

  /**
   * @brief Complete relocation after `ComputeStepAndNextVolume()`.
   * @param globalpoint Track endpoint in global coordinates.
   * @param globaldir Track direction in global coordinates.
   * @param[in,out] state Boundary state to descend through coincident volumes.
   */
  VECCORE_ATT_HOST_DEVICE
  static void RelocateToNextVolumeImpl(Vector3D<Precision> const &globalpoint, Vector3D<Precision> const &globaldir,
                                       vecgeom::NavigationState &state)
  {
    // if already outside, don't do anything
    if (state.IsOutside()) return;

    // Push the point inside the next volume.
    // A.G. This should not be needed now since LocatePointIn is boundary-aware
    Vector3D<Precision> pushed = globalpoint /* + kBoundaryPush * globaldir*/;

    // Calculate local point from global point.
    vecgeom::Transformation3D m;
    state.TopMatrix(m);
    Vector3D<Precision> localpoint = m.Transform(pushed);

    // passing the state to check in + the local point in the reference frame of the state
    LocatePointInNavState(localpoint, state, false, state.GetLastExited());

    if (state.Top() != nullptr) {
      while (state.Top()->IsAssembly()) {
        state.Pop();
      }
      VECGEOM_ASSERT(!state.Top()->GetLogicalVolume()->GetUnplacedVolume()->IsAssembly());
    }
  }

  /** @} */
};

/// Explicit name for the navigator backed by IndexedBVH.
using IndexedBVHNavigator = BVHNavigator;

} // namespace vecgeom

#endif // BVH_NAVIGATOR_H
