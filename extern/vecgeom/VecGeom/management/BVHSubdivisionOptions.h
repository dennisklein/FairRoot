//===-- VecGeom/management/BVHSubdivisionOptions.h ------------------*- C++ -*-===//
//
/// @file BVHSubdivisionOptions.h
/// @brief Shared host-side configuration selectors for BVH AABB subdivision.
//
//===----------------------------------------------------------------------===//

#ifndef VECGEOM_MANAGEMENT_BVHSUBDIVISIONOPTIONS_H_
#define VECGEOM_MANAGEMENT_BVHSUBDIVISIONOPTIONS_H_

#include "VecGeom/base/Global.h"

namespace vecgeom {

/// Largest supported number of uniform slabs for one daughter split.
///
/// Current automatic-mode defaults are at most 20. This higher structural
/// bound leaves room for expert studies while preventing accidental or
/// malformed inputs from requesting impractical allocations.
inline constexpr int kMaximumAABBSubdivisionM = 1024;

/// Build/navigation trade-off used by automatic daughter-AABB subdivision.
enum class AABBSubdivisionMode {
  kOff,        ///< Never split; use one AABB per daughter.
  kFast,       ///< Bounded clipped-extent search with low build cost.
  kBalanced,   ///< Moderate build/memory cost for stronger pruning.
  kNavigation, ///< Higher build/memory cost, optimized for navigation speed.
};

/// Configuration for automatic AABB subdivision.
///
/// The automatic modes are heuristic selectors for build/navigation trade-offs.
/// They never change the geometric candidate id: every returned slab remains a
/// conservative mother-frame AABB for the original placed volume. These knobs
/// are host-side BVH build inputs and must be set before BVHManager::Init().
struct AABBSubdividerConfig {
  /// Automatic search mode.
  AABBSubdivisionMode mode = AABBSubdivisionMode::kOff;

  /// Hard upper limit on M in navigation mode.
  int autoMaxM = 20;

  /// Maximum M in fast mode. This keeps build cost bounded.
  int fastMaxM = 4;

  /// Maximum M in balanced mode.
  int balancedMaxM = 8;

  /// Minimum unpenalized slab/original surface-area reduction required in fast mode.
  Precision fastMinReduction = Precision(0.15);

  /// Minimum transformed-AABB/local-extent surface-area ratio for fast mode.
  Precision fastMinAABBInflation = Precision(5.0);

  /// Minimum unpenalized slab/original surface-area reduction required in balanced mode.
  Precision balancedMinReduction = Precision(0.15);

  /// Minimum transformed-AABB/local-extent surface-area ratio tried by balanced mode.
  Precision balancedMinAABBInflation = Precision(1.5);

  /// Balanced-mode transformed/local AABB surface-area ratio above which all axes are tested.
  ///
  /// Candidate daughters with a lower ratio still use the cheaper one-axis
  /// search, so balanced mode spends the expensive all-axis search only on
  /// very inefficient transformed AABBs.
  Precision balancedFullSearchAABBInflation = Precision(3.0);

  /// Minimum unpenalized slab/original surface-area reduction required in navigation mode.
  Precision navigationMinReduction = Precision(0.02);

  /// Stop balanced/navigation M scans after this many non-improving M values.
  int nonImprovingStopCount = 2;

  /// Objective penalty per extra slab in fast mode, in units of SA(M=1).
  Precision fastBoxPenalty = Precision(0.030);

  /// Objective penalty per extra slab in balanced mode, in units of SA(M=1).
  Precision balancedBoxPenalty = Precision(0.020);

  /// Objective penalty per extra slab in navigation mode, in units of SA(M=1).
  Precision navigationBoxPenalty = Precision(0.002);

  /// Minimum estimated whole-mother BVH traversal-cost reduction required in fast mode; negative disables the gate.
  Precision fastMinCostReduction = Precision(0.35);

  /// Minimum estimated whole-mother BVH traversal-cost reduction required in balanced mode; negative disables the gate.
  Precision balancedMinCostReduction = Precision(0.20);

  /// Minimum estimated whole-mother BVH traversal-cost reduction required in navigation mode; negative disables the
  /// gate.
  Precision navigationMinCostReduction = Precision(0.20);

  /// Relative cost assigned to visiting one BVH node by the post-build gate.
  Precision nodeCost = Precision(1.0);

  /// Relative cost assigned to testing one primitive AABB by the post-build gate.
  Precision boxCost = Precision(1.0);

  /// Relative cost assigned to one daughter solid distance call by the post-build gate.
  Precision solidCost = Precision(5.0);

  /// Maximum mother daughter count eligible for subdivision in fast mode.
  ///
  /// Larger mothers keep identity-mapped BVHs. This bounds the temporary build
  /// products and avoids duplicating full-mother AABB/id arrays when only a
  /// sparse subset of daughters would benefit.
  int fastMaxDaughtersPerMother = 100000;

  /// Minimum fraction of a mother's daughters that must pass the fast prefilter.
  ///
  /// Requiring a representative fraction amortizes the mapped split-box
  /// representation over enough useful daughters. A value of 0 disables this
  /// fast-mode gate.
  Precision fastMinCandidateFractionPerMother = Precision(0.15);
};

/// Return a stable lower-case label for an AABB subdivision mode.
inline const char *AABBSubdivisionModeName(AABBSubdivisionMode mode)
{
  switch (mode) {
  case AABBSubdivisionMode::kOff:
    return "off";
  case AABBSubdivisionMode::kFast:
    return "fast";
  case AABBSubdivisionMode::kBalanced:
    return "balanced";
  case AABBSubdivisionMode::kNavigation:
    return "navigation";
  }
  return "unknown";
}

} // namespace vecgeom

#endif // VECGEOM_MANAGEMENT_BVHSUBDIVISIONOPTIONS_H_
