//===-- VecGeom/volumes/utilities/AABBSubdivider.h ------------------*- C++ -*-===//
//
/// @file   AABBSubdivider.h
/// @brief  Split a placed volume's bounding box into conservative clipped sub-boxes.
///
/// A single axis-aligned bounding box around an elongated, rotated or
/// diagonally-placed solid encloses empty space. AABBSubdivider replaces that
/// one loose box by a small collection of mother-frame slab boxes. Each slab is
/// computed from the solid's host-only CalculateClippedExtent() contract; no
/// sampling or Inside() probing is used here.
///
/// Main entry points:
///   Subdivide(placed, M)   -> exactly M slabs along the best axis.
///   SubdivideAuto(placed)  -> chooses M and axis from the configured
///                             build/navigation trade-off mode.
/// Both return boxes in the mother frame; M = 1 gives the plain transformed
/// Extent() AABB used by the current BVH baseline.
//
//===--------------------------------------------------------------------------===//

#pragma once

#include "VecGeom/base/Vector3D.h"
#include "VecGeom/base/Transformation3D.h"
#include "VecGeom/management/BVHSubdivisionOptions.h"
#include "VecGeom/volumes/PlacedVolume.h"
#include "VecGeom/volumes/LogicalVolume.h"
#include "VecGeom/volumes/UnplacedVolume.h"
#include "VecGeom/volumes/utilities/ClippedExtent.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

namespace vecgeom {
inline namespace cxx {

/// One conservative axis-aligned bounding box in the mother frame.
struct SlabBox {
  Vector3D<Precision> min{kInfLength, kInfLength, kInfLength};
  Vector3D<Precision> max{-kInfLength, -kInfLength, -kInfLength};
  int slabIndex{0}; ///< 0-based within the generated slabs.

  bool Valid() const { return min[0] <= max[0] && min[1] <= max[1] && min[2] <= max[2]; }
  Precision SurfaceArea() const
  {
    if (!Valid()) return Precision(0);
    auto d = max - min;
    return Precision(2) * (d[0] * d[1] + d[1] * d[2] + d[2] * d[0]);
  }
  Vector3D<Precision> HalfExtents() const { return (max - min) * Precision(0.5); }
  Vector3D<Precision> Center() const { return (max + min) * Precision(0.5); }
};

/// Diagnostic result returned by SubdivideAuto.
///
/// `slabs` is always conservative. When no accepted split satisfies the mode's
/// reduction and memory objective, the result contains the single transformed
/// local-extent AABB with `chosenM == 1`.
struct AABBSubdivisionResult {
  std::vector<SlabBox> slabs;       ///< Final conservative boxes in mother frame.
  int chosenM{1};                   ///< M that was selected.
  int chosenAxis{0};                ///< 0=X, 1=Y, 2=Z.
  Precision globalSA{0};            ///< Surface area of the M=1 global AABB.
  Precision finalSA{0};             ///< Total surface area of chosen slabs.
  bool usedSingleAxisSearch{false}; ///< Auto search fixed one axis instead of testing all axes.
  /// Per-step SA values indexed by M-1 (entry 0 = global SA at M=1).
  std::vector<Precision> saHistory;
};

class AABBSubdivider {
public:
  static AABBSubdivider &Instance()
  {
    static AABBSubdivider inst;
    return inst;
  }

  /// Subdivide @p placed into exactly @p M clipped-extent AABBs.
  ///
  /// @param placed Placed volume.
  /// @param M Number of uniform slabs along the selected axis. M=1 returns the
  ///          transformed Extent() AABB used by the BVH baseline. Values above
  ///          `kMaximumAABBSubdivisionM` return the unsplit fallback.
  /// @param cfg Tuning knobs for axis objective.
  /// @return Conservative slab boxes in mother frame.
  std::vector<SlabBox> Subdivide(VPlacedVolume const *placed, int M, AABBSubdividerConfig const &cfg = {}) const;

  /// @brief Subdivide @p placed into one previously selected uniform slab sequence.
  ///
  /// This entry point reconstructs a persisted split decision without testing
  /// other axes or values of M. Slab bounds are derived from the transformed,
  /// padded baseline AABB, so callers do not need to retain an offset or width.
  ///
  /// @param placed Placed volume.
  /// @param axisIdx Mother-frame split axis: 0=X, 1=Y, or 2=Z.
  /// @param M Number of uniform slabs. M=1 returns the unsplit baseline AABB;
  ///          values above `kMaximumAABBSubdivisionM` return no slabs.
  /// @return Conservative slab boxes in mother coordinates, or an empty vector
  ///         when no valid slab can be reconstructed.
  std::vector<SlabBox> SubdivideAlongAxis(VPlacedVolume const *placed, int axisIdx, int M) const;

  /// @brief Subdivide using an available padded mother-frame baseline AABB.
  ///
  /// @param placed Placed volume.
  /// @param axisIdx Mother-frame split axis: 0=X, 1=Y, or 2=Z.
  /// @param M Number of uniform slabs; values above
  ///          `kMaximumAABBSubdivisionM` return no slabs.
  /// @param globalBox Baseline AABB whose selected-axis bounds define the slab
  ///                  offset and spacing.
  /// @return Conservative slab boxes in mother coordinates, or an empty vector
  ///         when no valid slab can be reconstructed.
  /// @pre @p globalBox is a valid conservative padded AABB of @p placed in the
  ///      mother frame. M=1 returns this box unchanged.
  std::vector<SlabBox> SubdivideAlongAxis(VPlacedVolume const *placed, int axisIdx, int M,
                                          SlabBox const &globalBox) const;

  /// Automatically determine M and axis using the configured search mode.
  ///
  /// Fast mode chooses one split axis from the transformed AABB shape and scans
  /// a small M range on that axis with a strong memory/build penalty.
  /// Balanced mode is tiered: moderately inflated daughter boxes use the same
  /// one-axis search as fast mode, while highly inflated boxes test all axes.
  /// Navigation mode scans up to autoMaxM with a small penalty, favoring
  /// navigation speed and accepting more split boxes.
  ///
  /// @param placed Placed volume.
  /// @param cfg Tuning knobs and mode.
  /// @return Selected slab boxes plus diagnostics. chosenM=1 means no split.
  AABBSubdivisionResult SubdivideAuto(VPlacedVolume const *placed, AABBSubdividerConfig const &cfg = {}) const;

private:
  /// Return the unsplit transformed local-extent AABB used as the baseline.
  static SlabBox ComputeGlobalAABB(VUnplacedVolume const *unplaced, Transformation3D const *tr);

  /// Compute one conservative mother-frame slab box for a voxel limit.
  ///
  /// If a solid returns inconsistent existence across axes, the helper falls
  /// back to clipping the unsplit AABB to avoid under-covering the daughter.
  static SlabBox ComputeOneSlabAABB(VUnplacedVolume const *unplaced, Transformation3D const *tr, int axisIdx,
                                    Precision lo, Precision hi, SlabBox const &globalBox, int slabIndex);

  /// Build all valid uniform slabs along one selected mother-frame axis.
  static std::vector<SlabBox> ComputeSlabs(VUnplacedVolume const *unplaced, Transformation3D const *tr, int axisIdx,
                                           int M, SlabBox const &globalBox);

  /// Try all axes for a fixed M and return the lowest penalized objective.
  static std::vector<SlabBox> BestAxisSlabs(VUnplacedVolume const *unplaced, Transformation3D const *tr, int M,
                                            SlabBox const &globalBox, AABBSubdividerConfig const &cfg, int &outBestAxis,
                                            Precision &outSA);

  static Precision TotalSurfaceArea(std::vector<SlabBox> const &slabs);
  static Precision ExtentSurfaceArea(VUnplacedVolume const *unplaced);
  static Precision AABBInflationRatio(VUnplacedVolume const *unplaced, SlabBox const &globalBox);
  static Precision CandidateObjective(Precision totalSA, std::size_t nSlabs, Precision globalSA,
                                      AABBSubdividerConfig const &cfg);
  static Precision ModeMinReduction(AABBSubdividerConfig const &cfg);
  static int ModeMaxM(AABBSubdividerConfig const &cfg);
  static int LongestAxis(SlabBox const &box);
  static int EstimateFastMaxM(SlabBox const &box, int maxM);
};

inline SlabBox AABBSubdivider::ComputeGlobalAABB(VUnplacedVolume const *unplaced, Transformation3D const *tr)
{
  Vector3D<Precision> localMin, localMax;
  unplaced->Extent(localMin, localMax);
  const ClippedBox box = TransformLocalAABBToMother(localMin, localMax, tr);

  SlabBox result;
  result.min       = box.min;
  result.max       = box.max;
  result.slabIndex = 0;
  return result;
}

inline SlabBox AABBSubdivider::ComputeOneSlabAABB(VUnplacedVolume const *unplaced, Transformation3D const *tr,
                                                  int axisIdx, Precision lo, Precision hi, SlabBox const &globalBox,
                                                  int slabIndex)
{
  SlabBox result;
  result.slabIndex = slabIndex;

  VoxelLimits limits;
  limits.AddLimit(AxisFromIndex(axisIdx), lo, hi);

  bool axisResult[3]    = {false, false, false};
  Precision minValue[3] = {0., 0., 0.};
  Precision maxValue[3] = {0., 0., 0.};
  int nExisting         = 0;
  for (int axis = 0; axis < 3; ++axis) {
    axisResult[axis] =
        unplaced->CalculateClippedExtent(AxisFromIndex(axis), limits, *tr, minValue[axis], maxValue[axis]);
    if (axisResult[axis]) ++nExisting;
  }

  if (nExisting == 0) return result;

  if (nExisting == 3) {
    result.min = Vector3D<Precision>(minValue[0], minValue[1], minValue[2]);
    result.max = Vector3D<Precision>(maxValue[0], maxValue[1], maxValue[2]);
    return result;
  }

  // Projection existence should be identical for all axes. If a solid override
  // violates that contract, keep the slab conservative by clipping the global
  // AABB instead of trusting the partial answer.
  ClippedBox fallback(globalBox.min, globalBox.max);
  if (!IntersectBoxWithLimits(fallback, limits)) return result;
  result.min = fallback.min;
  result.max = fallback.max;
  return result;
}

inline std::vector<SlabBox> AABBSubdivider::ComputeSlabs(VUnplacedVolume const *unplaced, Transformation3D const *tr,
                                                         int axisIdx, int M, SlabBox const &globalBox)
{
  const Precision lo0 = globalBox.min[axisIdx];
  const Precision hi0 = globalBox.max[axisIdx];
  if (!(hi0 > lo0) || M < 1 || M > kMaximumAABBSubdivisionM) return {};

  const Precision slabW = (hi0 - lo0) / Precision(M);
  std::vector<SlabBox> slabs;
  slabs.reserve(M);
  for (int k = 0; k < M; ++k) {
    auto slab = ComputeOneSlabAABB(unplaced, tr, axisIdx, lo0 + k * slabW, lo0 + (k + 1) * slabW, globalBox, k);
    if (slab.Valid()) slabs.push_back(slab);
  }
  return slabs;
}

inline Precision AABBSubdivider::TotalSurfaceArea(std::vector<SlabBox> const &slabs)
{
  Precision total = 0.;
  for (auto const &slab : slabs)
    total += slab.SurfaceArea();
  return total;
}

inline Precision AABBSubdivider::ExtentSurfaceArea(VUnplacedVolume const *unplaced)
{
  Vector3D<Precision> localMin, localMax;
  unplaced->Extent(localMin, localMax);
  SlabBox localBox;
  localBox.min = localMin;
  localBox.max = localMax;
  return localBox.SurfaceArea();
}

inline Precision AABBSubdivider::AABBInflationRatio(VUnplacedVolume const *unplaced, SlabBox const &globalBox)
{
  const Precision localSA = ExtentSurfaceArea(unplaced);
  return localSA > 0. ? globalBox.SurfaceArea() / localSA : Precision(0);
}

inline Precision AABBSubdivider::CandidateObjective(Precision totalSA, std::size_t nSlabs, Precision globalSA,
                                                    AABBSubdividerConfig const &cfg)
{
  if (globalSA <= 0.) return std::numeric_limits<Precision>::max();
  Precision penalty = cfg.fastBoxPenalty;
  if (cfg.mode == AABBSubdivisionMode::kBalanced) penalty = cfg.balancedBoxPenalty;
  if (cfg.mode == AABBSubdivisionMode::kNavigation) penalty = cfg.navigationBoxPenalty;
  return totalSA / globalSA + penalty * Precision(nSlabs > 0 ? nSlabs - 1 : 0);
}

inline Precision AABBSubdivider::ModeMinReduction(AABBSubdividerConfig const &cfg)
{
  if (cfg.mode == AABBSubdivisionMode::kNavigation) return cfg.navigationMinReduction;
  if (cfg.mode == AABBSubdivisionMode::kBalanced) return cfg.balancedMinReduction;
  return cfg.fastMinReduction;
}

inline int AABBSubdivider::ModeMaxM(AABBSubdividerConfig const &cfg)
{
  if (cfg.mode == AABBSubdivisionMode::kNavigation) return std::clamp(cfg.autoMaxM, 1, kMaximumAABBSubdivisionM);
  if (cfg.mode == AABBSubdivisionMode::kBalanced) return std::clamp(cfg.balancedMaxM, 1, kMaximumAABBSubdivisionM);
  if (cfg.mode == AABBSubdivisionMode::kFast) return std::clamp(cfg.fastMaxM, 1, kMaximumAABBSubdivisionM);
  return 1;
}

inline int AABBSubdivider::LongestAxis(SlabBox const &box)
{
  const auto extent = box.max - box.min;
  int axis          = 0;
  if (extent[1] > extent[axis]) axis = 1;
  if (extent[2] > extent[axis]) axis = 2;
  return axis;
}

inline int AABBSubdivider::EstimateFastMaxM(SlabBox const &box, int maxM)
{
  if (maxM <= 1) return 1;

  // Fast mode uses only the transformed AABB aspect ratio. This intentionally
  // avoids trying several expensive clipped-extent scenarios before deciding
  // whether a production-style split is worth building.
  Precision extent[3] = {box.max[0] - box.min[0], box.max[1] - box.min[1], box.max[2] - box.min[2]};
  std::sort(extent, extent + 3, [](Precision a, Precision b) { return a > b; });
  if (!(extent[1] > 0.)) return 1;

  const int aspectM = static_cast<int>(std::ceil(extent[0] / extent[1]));
  return std::max(2, std::min(maxM, aspectM));
}

inline std::vector<SlabBox> AABBSubdivider::BestAxisSlabs(VUnplacedVolume const *unplaced, Transformation3D const *tr,
                                                          int M, SlabBox const &globalBox,
                                                          AABBSubdividerConfig const &cfg, int &outBestAxis,
                                                          Precision &outSA)
{
  std::vector<SlabBox> bestSlabs;
  Precision bestObjective = std::numeric_limits<Precision>::max();
  outBestAxis             = 0;
  outSA                   = std::numeric_limits<Precision>::max();

  const Precision globalSA = globalBox.SurfaceArea();
  for (int axis = 0; axis < 3; ++axis) {
    if (globalBox.max[axis] - globalBox.min[axis] < Precision(1e-10)) continue;

    auto slabs                  = ComputeSlabs(unplaced, tr, axis, M, globalBox);
    const Precision totalSA     = TotalSurfaceArea(slabs);
    const Precision objective   = CandidateObjective(totalSA, slabs.size(), globalSA, cfg);
    const bool betterObjective  = objective < bestObjective;
    const bool tiedButTighterSA = objective == bestObjective && totalSA < outSA;
    if (betterObjective || tiedButTighterSA) {
      bestObjective = objective;
      outBestAxis   = axis;
      outSA         = totalSA;
      bestSlabs     = std::move(slabs);
    }
  }
  return bestSlabs;
}

inline std::vector<SlabBox> AABBSubdivider::Subdivide(VPlacedVolume const *placed, int M,
                                                      AABBSubdividerConfig const &cfg) const
{
  assert(placed && "AABBSubdivider::Subdivide: null placed volume");
  assert(M >= 1 && "AABBSubdivider::Subdivide: M must be >= 1");

  VUnplacedVolume const *unplaced = placed->GetLogicalVolume()->GetUnplacedVolume();
  Transformation3D const *tr      = placed->GetTransformation();
  SlabBox globalBox               = ComputeGlobalAABB(unplaced, tr);
  if (M == 1) return {globalBox};

  int bestAxis      = 0;
  Precision totalSA = 0.;
  auto slabs        = BestAxisSlabs(unplaced, tr, M, globalBox, cfg, bestAxis, totalSA);
  return slabs.empty() ? std::vector<SlabBox>{globalBox} : slabs;
}

inline std::vector<SlabBox> AABBSubdivider::SubdivideAlongAxis(VPlacedVolume const *placed, int axisIdx, int M) const
{
  assert(placed && "AABBSubdivider::SubdivideAlongAxis: null placed volume");
  assert(axisIdx >= 0 && axisIdx < 3 && "AABBSubdivider::SubdivideAlongAxis: invalid axis");
  assert(M >= 1 && "AABBSubdivider::SubdivideAlongAxis: M must be >= 1");

  VUnplacedVolume const *unplaced = placed->GetLogicalVolume()->GetUnplacedVolume();
  Transformation3D const *tr      = placed->GetTransformation();
  SlabBox globalBox               = ComputeGlobalAABB(unplaced, tr);

  return SubdivideAlongAxis(placed, axisIdx, M, globalBox);
}

inline std::vector<SlabBox> AABBSubdivider::SubdivideAlongAxis(VPlacedVolume const *placed, int axisIdx, int M,
                                                               SlabBox const &globalBox) const
{
  assert(placed && "AABBSubdivider::SubdivideAlongAxis: null placed volume");
  assert(axisIdx >= 0 && axisIdx < 3 && "AABBSubdivider::SubdivideAlongAxis: invalid axis");
  assert(M >= 1 && "AABBSubdivider::SubdivideAlongAxis: M must be >= 1");
  assert(globalBox.Valid() && "AABBSubdivider::SubdivideAlongAxis: invalid baseline AABB");

  if (M == 1) return {globalBox};

  VUnplacedVolume const *unplaced = placed->GetLogicalVolume()->GetUnplacedVolume();
  Transformation3D const *tr      = placed->GetTransformation();
  return ComputeSlabs(unplaced, tr, axisIdx, M, globalBox);
}

inline AABBSubdivisionResult AABBSubdivider::SubdivideAuto(VPlacedVolume const *placed,
                                                           AABBSubdividerConfig const &cfg) const
{
  assert(placed && "AABBSubdivider::SubdivideAuto: null placed volume");

  VUnplacedVolume const *unplaced = placed->GetLogicalVolume()->GetUnplacedVolume();
  Transformation3D const *tr      = placed->GetTransformation();
  SlabBox globalBox               = ComputeGlobalAABB(unplaced, tr);

  AABBSubdivisionResult res;
  res.slabs      = {globalBox};
  res.chosenM    = 1;
  res.chosenAxis = 0;
  res.globalSA   = globalBox.SurfaceArea();
  res.finalSA    = res.globalSA;
  res.saHistory.push_back(res.globalSA);

  if (res.globalSA < std::numeric_limits<Precision>::epsilon() || cfg.mode == AABBSubdivisionMode::kOff) return res;

  std::vector<SlabBox> bestSlabs = res.slabs;
  Precision bestSA               = res.globalSA;
  Precision bestObjective        = CandidateObjective(bestSA, bestSlabs.size(), res.globalSA, cfg);
  int bestAxis                   = 0;
  int nonImproving               = 0;
  const int maxM                 = ModeMaxM(cfg);
  const bool fastMode            = cfg.mode == AABBSubdivisionMode::kFast;
  // Balanced mode is tiered: moderately inflated daughters use the same cheap
  // one-axis search as fast mode, while only the loosest boxes pay all-axis
  // search cost.
  const bool balancedSingleAxis = cfg.mode == AABBSubdivisionMode::kBalanced &&
                                  AABBInflationRatio(unplaced, globalBox) <= cfg.balancedFullSearchAABBInflation;
  const bool singleAxisMode     = fastMode || balancedSingleAxis;
  res.usedSingleAxisSearch      = singleAxisMode;
  const int fastAxis            = singleAxisMode ? LongestAxis(globalBox) : 0;
  const int scanMaxM            = singleAxisMode ? EstimateFastMaxM(globalBox, maxM) : maxM;

  for (int M = 2; M <= scanMaxM; ++M) {
    int candidateAxis     = 0;
    Precision candidateSA = std::numeric_limits<Precision>::max();
    // Single-axis tiers commit to one cheap axis estimate. Full balanced and
    // navigation tiers pay the higher build cost of testing all axes for each M.
    auto candidateSlabs = singleAxisMode ? ComputeSlabs(unplaced, tr, fastAxis, M, globalBox)
                                         : BestAxisSlabs(unplaced, tr, M, globalBox, cfg, candidateAxis, candidateSA);
    if (singleAxisMode) {
      candidateAxis = fastAxis;
      candidateSA   = TotalSurfaceArea(candidateSlabs);
    }
    res.saHistory.push_back(candidateSA);
    if (candidateSlabs.empty()) {
      ++nonImproving;
    } else {
      const Precision objective = CandidateObjective(candidateSA, candidateSlabs.size(), res.globalSA, cfg);
      if (objective < bestObjective) {
        bestSlabs     = std::move(candidateSlabs);
        bestSA        = candidateSA;
        bestObjective = objective;
        bestAxis      = candidateAxis;
        res.chosenM   = M;
        nonImproving  = 0;
      } else {
        ++nonImproving;
      }
    }

    if (!fastMode && nonImproving >= cfg.nonImprovingStopCount) break;
  }

  const Precision rawReduction = (res.globalSA - bestSA) / res.globalSA;
  if (res.chosenM <= 1 || rawReduction < ModeMinReduction(cfg)) return res;

  res.slabs      = std::move(bestSlabs);
  res.chosenAxis = bestAxis;
  res.finalSA    = bestSA;
  return res;
}

} // namespace cxx
} // namespace vecgeom
