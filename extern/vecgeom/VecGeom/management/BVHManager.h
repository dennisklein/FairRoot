/// \file BVHManager.h
/// \author Guilherme Amadio

#ifndef VECGEOM_MANAGEMENT_BVHMANAGER_H_
#define VECGEOM_MANAGEMENT_BVHMANAGER_H_

#include "VecGeom/base/Assert.h"
#include "VecGeom/base/Config.h"
#include "VecGeom/base/Cuda.h"
#include "VecGeom/base/IndexedBVH.h"
#include "VecGeom/management/BVHSubdivisionOptions.h"
#include "VecGeom/volumes/LogicalVolume.h"

#include <cstddef>
#include <string>

#ifdef VECGEOM_ENABLE_CUDA
#include "VecGeom/backend/cuda/Interface.h"
#endif

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {
template <typename Real_t>
inline std::vector<IndexedBVH<Real_t> *> hBVH;
#ifdef VECGEOM_ENABLE_CUDA
template <typename Real_t>
inline __device__ IndexedBVH<Real_t> *dBVH;
#endif

// Macro allowing downstream codes to use GetDeviceBVH
#define VECGEOM_BVHMANAGER_DEVICE

/// Host-side diagnostics collected during the most recent BVHManager build.
///
/// Timings reported from parallel build phases distinguish wall time from
/// accumulated worker time. Accumulated worker time may exceed initElapsed when
/// several workers execute clipped-extent preparation or BVH construction
/// concurrently.
struct BVHBuildStatistics {
  double initElapsed                   = 0.; ///< Total BVHManager::Init() wall time.
  double boxPreparationElapsed         = 0.; ///< Wall time spent computing per-volume BVH box products.
  double candidateSelectionWallElapsed = 0.; ///< Wall time spent running the parallel daughter prefilter.
  double candidateSelectionElapsed     = 0.; ///< Accumulated worker time spent in split prefilter checks.
  double candidatePreparationElapsed   = 0.; ///< Serial time spent reducing prefilters into candidate jobs.
  double subdivisionWallElapsed        = 0.; ///< Wall time spent running parallel clipped-extent subdivision.
  double subdivisionElapsed           = 0.; ///< Accumulated worker time spent building clipped-extent split candidates.
  double boxAssemblyElapsed           = 0.; ///< Serial time spent assembling daughter products by mother.
  double bvhConstructionElapsed       = 0.; ///< Wall time spent constructing IndexedBVHs.
  double bvhConstructionWorkerElapsed = 0.; ///< Accumulated worker time constructing IndexedBVHs.
  double deviceInitElapsed            = 0.; ///< Time spent copying IndexedBVHs to device storage.
  double indexedBVHConstructionElapsed       = 0.; ///< IndexedBVH construction wall time.
  double indexedBVHConstructionWorkerElapsed = 0.; ///< Accumulated IndexedBVH worker time.
  double indexedDeviceInitElapsed            = 0.; ///< IndexedBVH device-copy time.
  double recipeValidationElapsed      = 0.; ///< Wall time spent validating and matching persisted split decisions.
  double recipeRebuildWallElapsed     = 0.; ///< Wall time spent rebuilding the selected persisted slabs.
  double recipeRebuildElapsed         = 0.; ///< Accumulated worker time spent rebuilding persisted slabs.
  double recipeAutoWriteElapsed       = 0.; ///< Wall time spent preparing and publishing an absent environment recipe.
  std::size_t logicalVolumes          = 0;  ///< Logical volumes visited by Init().
  std::size_t logicalVolumesWithBoxes = 0;  ///< Logical volumes with at least one daughter AABB.
  std::size_t splitLogicalVolumes     = 0;  ///< Logical volumes where at least one daughter was split.
  std::size_t originalBoxes           = 0;  ///< Number of original daughter AABBs.
  std::size_t finalBoxes              = 0;  ///< Number of BVH primitive boxes after optional subdivision.
  std::size_t splitCandidates         = 0;  ///< Daughters passing the split prefilter.
  std::size_t balancedSingleAxisCandidates = 0;  ///< Balanced daughters using one-axis split search.
  std::size_t balancedMultiAxisCandidates  = 0;  ///< Balanced daughters using all-axis split search.
  std::size_t acceptedSplits               = 0;  ///< Daughters where subdivision was accepted.
  std::size_t splitBoxes                   = 0;  ///< BVH boxes produced by accepted subdivisions.
  std::size_t extraBoxes                   = 0;  ///< Added boxes, i.e. finalBoxes - originalBoxes for split daughters.
  std::size_t costEvaluatedLogicalVolumes  = 0;  ///< Split mothers evaluated by the post-build BVH cost gate.
  std::size_t costRejectedLogicalVolumes   = 0;  ///< Split mothers reverted to their original BVH by the cost gate.
  std::size_t costRejectedSplits           = 0;  ///< Daughter split decisions removed by the cost gate.
  std::size_t costRejectedExtraBoxes       = 0;  ///< Added primitive boxes removed by the cost gate.
  std::size_t fastMotherSkippedDaughters   = 0;  ///< Daughters skipped by the fast-mode mother policy.
  std::size_t fastInflationSamples         = 0;  ///< Daughters with a valid original/local AABB SA ratio.
  double fastInflationSum                  = 0.; ///< Sum of original/local AABB surface-area ratios.
  double fastInflationMax                  = 0.; ///< Maximum original/local AABB surface-area ratio.
  std::size_t fastInflationAbove1p5        = 0;  ///< Daughters with original/local ratio > 1.5.
  std::size_t fastInflationAbove2          = 0;  ///< Daughters with original/local ratio > 2.
  std::size_t fastInflationAbove3          = 0;  ///< Daughters with original/local ratio > 3.
  std::size_t fastInflationAbove5          = 0;  ///< Daughters with original/local ratio > 5.
  std::size_t fastInflationAbove8          = 0;  ///< Daughters with original/local ratio > 8.
  std::size_t fastInflationAbove13         = 0;  ///< Daughters with original/local ratio > 13.
  std::size_t recipeEntriesLoaded          = 0;  ///< Logical-volume entries parsed from the active recipe file.
  std::size_t recipeEntriesUsed            = 0;  ///< Fresh entries used without live heuristic evaluation.
  std::size_t recipeEntriesRejected        = 0;  ///< Incompatible, stale, malformed-product, or orphan entries.
  std::size_t recipeEntriesMissing         = 0;  ///< Current mothers absent from an otherwise compatible recipe file.
  std::size_t recipeEntriesLiveBuilt       = 0;  ///< Current mothers sent through the live fallback path.
  std::size_t recipeSplitDecisionsUsed     = 0;  ///< Persisted fixed-axis/M daughter decisions used.
  std::size_t recipeSplitBoxesBuilt        = 0;  ///< Final boxes rebuilt from persisted decisions.
  std::size_t recipeDecisionBytesRetained  = 0;  ///< Approximate host bytes retained for later recipe export.
  std::size_t recipeAutoWriteCreated       = 0;  ///< Missing environment recipe files created after the build.
  std::size_t recipeAutoWriteSkipped       = 0;  ///< Recipe files published first by another initializer.
  std::size_t recipeAutoWriteFailed        = 0;  ///< Missing environment recipe files that could not be published.
  std::size_t hostAllocatedBytes           = 0;  ///< Host arrays owned by IndexedBVHs.
  std::size_t deviceAllocatedBytes         = 0;  ///< Packed IndexedBVH device arena, including object slots.
  std::size_t deviceAllocationCount        = 0;  ///< Persistent device allocations for the active BVH table.
  std::size_t indexedHostAllocatedBytes    = 0;  ///< Bytes in IndexedBVH host arrays.
  std::size_t indexedDeviceAllocatedBytes  = 0;  ///< Bytes in the packed IndexedBVH device arena.

  std::size_t distanceCostRejectedLogicalVolumes = 0; ///< Split mothers reverted by the ray/Distance cost gate.
};

/// Range callback used by BVH build executors.
///
/// @param begin First build-job index in the range.
/// @param end One-past-last build-job index in the range.
/// @param context Opaque VecGeom-owned build context passed to ParallelFor().
using BVHBuildRangeTask = void (*)(std::size_t begin, std::size_t end, void *context);

/// Synchronous scheduling interface for host BVH construction.
///
/// Applications that already own geometry-initialization parallelism can provide
/// an adapter that dispatches VecGeom build ranges through their scheduler. The
/// executor must cover every index in [begin, end) by invoking @p task on one
/// or more subranges, must not retain @p task or @p context after returning,
/// and must propagate task failures to the caller. VecGeom owns the meaning of
/// the indices and the work performed by the callback. One Init() call may
/// issue multiple synchronous ranges for distinct build phases, such as cheap
/// candidate prefiltering, clipped-extent subdivision, and independent
/// logical-volume BVH construction.
class BVHBuildExecutor {
public:
  virtual ~BVHBuildExecutor() = default;

  /// Schedule VecGeom BVH build jobs.
  ///
  /// @param begin First build-job index.
  /// @param end One-past-last build-job index.
  /// @param minGrainSize Suggested minimum scheduler chunk size.
  /// @param task VecGeom callback to execute for each assigned subrange.
  /// @param context Opaque VecGeom-owned callback context.
  virtual void ParallelFor(std::size_t begin, std::size_t end, std::size_t minGrainSize, BVHBuildRangeTask task,
                           void *context) = 0;
};

/**
 * @brief The @c BVHManager class is a singleton class to manage the association between
 * logical volumes and their bounding volume hierarchies, using the logical volumes' ids.
 */
class BVHManager {

// to avoid changing the user code via templates, the precision is chosen at compile-time in the BVHManager and
// BVHNavigator
#ifdef VECGEOM_BVH_SINGLE
  using Real_t = float;
#else
  using Real_t = double;
#endif

public:
  BVHManager() = delete;

  /// Sentinel selecting the standard library's hardware-concurrency count.
  ///
  /// The numeric value is intentionally zero so APIs and command-line tools can
  /// accept a compact sentinel while call sites use the self-documenting name.
  static constexpr unsigned kCPUCount = 0;

  /// Select the global daughter-AABB subdivision mode used by Init().
  ///
  /// @param mode Build/navigation trade-off to use for future BVH construction.
  ///
  /// The setting must be applied before the geometry is closed or before Init()
  /// is called manually. Existing BVHs are not rebuilt by changing this setting.
  static void SetAABBSubdivisionMode(AABBSubdivisionMode mode);

  /// Return the currently configured global daughter-AABB subdivision mode.
  ///
  /// @return Selected mode. On host builds the value may be initialized from
  ///         VECGEOM_BVH_AABB_SUBDIVISION unless explicitly set through the API.
  static AABBSubdivisionMode GetAABBSubdivisionMode();

  /// Override the maximum number of slabs considered by automatic subdivision.
  ///
  /// @param maxM Positive maximum M to use for every mode, capped at
  ///             `kMaximumAABBSubdivisionM`, or <=0 to restore each mode's default.
  ///
  /// The setting must be applied before the geometry is closed or before Init()
  /// is called manually. Existing BVHs are not rebuilt by changing this setting.
  static void SetAABBSubdivisionMaxM(int maxM);

  /// Return the configured maximum M override, or 0 when mode defaults are used.
  ///
  /// @return Positive M override, or 0 when the selected mode uses its default.
  static int GetAABBSubdivisionMaxM();

  /// Replace the full host-side subdivision configuration used by Init().
  ///
  /// @param config Complete mode, M limits, split prefilter and mother-policy
  ///               triggers, reduction thresholds, and memory objective
  ///               penalties for future BVH construction.
  ///
  /// This call has API precedence over environment defaults. It also clears the
  /// legacy uniform max-M override, so the per-mode maxima in @p config are used
  /// unless SetAABBSubdivisionMaxM() is called afterwards.
  static void SetAABBSubdivisionConfig(AABBSubdividerConfig const &config);

  /// Return the effective host-side subdivision configuration.
  ///
  /// @return Current configuration after applying environment defaults and the
  ///         legacy uniform max-M override. Device compilation returns an `off`
  ///         configuration because subdivision is a host BVH build feature.
  static AABBSubdividerConfig GetAABBSubdivisionConfig();

  /// Return diagnostics from the most recent host BVH build.
  ///
  /// @return Snapshot of counters and timings collected by the last Init() call.
  static BVHBuildStatistics GetBuildStatistics();

  /// @brief Load a persisted daughter-AABB subdivision recipe for the next BVH build.
  ///
  /// @param fileName Recipe file to parse.
  /// @param[out] errorMessage Optional parse or I/O diagnostic; cleared on success.
  /// @return True when a complete supported recipe was loaded.
  ///
  /// Geometry freshness and storage-precision compatibility are checked by
  /// Init(), after baseline daughter AABBs are available. A valid decision takes
  /// precedence over the selected live heuristic and reconstructs only its
  /// selected fixed-axis, fixed-M slabs. A failed load clears any previously
  /// loaded recipe; stale or missing mothers use the live mode. The parsed
  /// payload is consumed and released by the next Init().
  static bool ReadAABBSubdivisionRecipe(std::string const &fileName, std::string *errorMessage = nullptr);

  /// @brief Persist split decisions from the most recent successful BVH build.
  ///
  /// @param fileName Recipe file to create or replace.
  /// @param[out] errorMessage Optional validation or I/O diagnostic; cleared on success.
  /// @return True when the complete recipe was written.
  ///
  /// The current geometry must be the geometry used by the most recent Init().
  /// Only accepted daughter splits are stored, as daughter id, selected
  /// mother-frame axis, and uniform slab count M. Unlisted daughters are
  /// intentional unsplit decisions. Loading reconstructs each listed split
  /// once, without testing alternative axes or values of M.
  static bool WriteAABBSubdivisionRecipe(std::string const &fileName, std::string *errorMessage = nullptr);

  /// @brief Clear the loaded recipe so future builds use only the configured live mode.
  static void ClearAABBSubdivisionRecipe();

  /// Install a non-owned synchronous executor for host BVH construction.
  ///
  /// @param executor Executor to use for future Init() calls, or nullptr to
  ///                 restore VecGeom's built-in executor.
  ///
  /// VecGeom stores this non-owned pointer until it is replaced or cleared, so
  /// the executor object must remain alive for every Init() call that can
  /// observe it. The built-in executor uses the configured thread cap and the
  /// same synchronous callback contract.
  static void SetBuildExecutor(BVHBuildExecutor *executor);

  /// Select the thread cap for VecGeom's built-in BVH build executor.
  ///
  /// @param nThreads `kCPUCount` for standard-library hardware concurrency, 1
  ///                 for forced serial execution, or N>1 to cap workers.
  ///
  /// The setting applies to future Init() calls. It is ignored while an
  /// external BVHBuildExecutor is installed through SetBuildExecutor().
  static void SetBuildThreadCount(unsigned nThreads);

  /// Return the configured built-in BVH build thread cap.
  ///
  /// @return Requested cap: `kCPUCount`, 1, or an explicit N>1 value. This is
  ///         not the resolved number of workers for a specific geometry.
  static unsigned GetBuildThreadCount();

  /**
   * Initializes the bounding volume hierarchies for all logical volumes in the geometry.
   * Since it uses the ABBoxManager to fetch the pre-computed bounding boxes for each logical volume,
   * it must be called after the bounding boxes have already been computed. The depth is not specified,
   * to allow the BVH class to choose the depth dynamically based on the number of children of each
   * logical volume.
   *
   * Call this explicitly after GeoManager::CloseGeometry(), which computes the
   * daughter AABBs consumed here. The BVHManager assumes all volumes have an
   * associated IndexedBVH, but only BVHs for volumes whose navigator is set to
   * BVHNavigator are actually accessed at runtime.
   *
   * Construction is transactional: a failed rebuild leaves the previously
   * published host BVHs available, while a successful rebuild replaces and
   * releases them only after every new hierarchy is complete.
   */
  static void Init();

  /** Initializes bounding volume hierarchies on the GPU. */
  static cuda::IndexedBVH<Real_t> const *DeviceInit();

  /** Access the device BVH pointer if CUDA is enabled. **/
  static cuda::IndexedBVH<Real_t> const *GetDeviceBVH();

  /** Compatibility name for the IndexedBVH device array. */
  static cuda::IndexedBVH<Real_t> const *GetDeviceIndexedBVH() { return GetDeviceBVH(); }

  VECCORE_ATT_HOST_DEVICE
  static IndexedBVH<Real_t> const *GetBVH(int id)
  {
#ifdef VECCORE_CUDA_DEVICE_COMPILATION
    return &cuda::dBVH<Real_t>[id];
#else
    // Init() fills the table; before it runs there is nothing to index.
    VECGEOM_ASSERT(id >= 0 && static_cast<size_t>(id) < hBVH<Real_t>.size());
    return hBVH<Real_t>[id];
#endif
  }

  VECCORE_ATT_HOST_DEVICE
  static IndexedBVH<Real_t> const *GetBVH(LogicalVolume const *v) { return GetBVH(v->id()); }

  /** Compatibility name for the active IndexedBVH. */
  VECCORE_ATT_HOST_DEVICE
  static IndexedBVH<Real_t> const *GetIndexedBVH(int id) { return GetBVH(id); }

  /** Compatibility overload accepting a logical-volume pointer. */
  VECCORE_ATT_HOST_DEVICE
  static IndexedBVH<Real_t> const *GetIndexedBVH(LogicalVolume const *v) { return GetBVH(v); }
};

} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif
