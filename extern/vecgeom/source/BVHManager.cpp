/// \file BVHManager.cpp
/// \author Guilherme Amadio

#include "VecGeom/management/BVHManager.h"
#include "VecGeom/management/ABBoxManager.h"
#ifdef VECGEOM_CUDA_INTERFACE
#include "VecGeom/management/CudaManager.h"
#endif
#include "VecGeom/management/Environment.h"
#include "VecGeom/management/GeoManager.h"
#include "VecGeom/management/Logger.h"
#ifndef VECCORE_CUDA
#include "BVHSubdivisionRecipe.h"
#include "VecGeom/base/Version.h"
#include "VecGeom/volumes/utilities/AABBSubdivider.h"
#endif

#include <algorithm>
#include <cassert>
#include <atomic>
#include <chrono>
#include <cctype>
#include <cstring>
#include <cstdlib>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <limits>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

namespace vecgeom {
#ifdef VECCORE_CUDA
inline
#endif
    namespace cuda {

template <typename Real_t>
IndexedBVH<Real_t> *InitializeDeviceIndexedBVHBuffer(void *arena, size_t n, size_t objectStride, size_t layoutSize);
template <typename Real_t>
IndexedBVH<Real_t> *GetDeviceIndexedBVHBuffer();

} // namespace cuda

inline namespace VECGEOM_IMPL_NAMESPACE {
#ifndef VECCORE_CUDA
namespace {

constexpr char const *kSubdivisionModeEnv   = "VECGEOM_BVH_AABB_SUBDIVISION";
constexpr char const *kSubdivisionMaxMEnv   = "VECGEOM_BVH_AABB_SUBDIVISION_MAX_M";
constexpr char const *kSubdivisionRecipeEnv = "VECGEOM_BVH_AABB_SUBDIVISION_RECIPE";

#ifdef VECGEOM_BVH_SINGLE
using BVHBuildReal_t = float;
#else
using BVHBuildReal_t = double;
#endif

using IndexedBVHOwnerTable = std::vector<std::unique_ptr<IndexedBVH<BVHBuildReal_t>>>;

AABBSubdividerConfig gSubdivisionConfig;
int gSubdivisionMaxM               = 0;
bool gSubdivisionModeExplicit      = false;
bool gSubdivisionMaxMExplicit      = false;
bool gSubdivisionEnvironmentLoaded = false;
bool gSubdivisionRecipeExplicit    = false;
std::string gSubdivisionRecipeEnvironmentFile;
std::string gSubdivisionRecipeCreateFile;
BVHBuildStatistics gBuildStats;
BVHBuildExecutor *gBuildExecutor = nullptr;
unsigned gBuildThreadCount       = BVHManager::kCPUCount;
std::unique_ptr<bvhrecipe::Recipe> gLoadedRecipe;
VPlacedVolume const *gLastBuildWorld = nullptr;
bool gHasSuccessfulBuild             = false;
std::vector<LogicalVolume const *> gLastBuildLogicalVolumes;
// hBVH remains the raw lookup table used by navigation. This parallel owner table
// makes successful rebuilds release the previous hierarchy and also cleans up
// the final hierarchy during process shutdown.
IndexedBVHOwnerTable gOwnedIndexedBVHs;
// Retain only accepted decisions so post-Init export does not keep a second
// copy of final slab coordinates or empty per-logical-volume vectors.
std::unordered_map<LogicalVolume const *, std::vector<bvhrecipe::Daughter>> gLastBuildSplitDecisions;

using Clock = std::chrono::steady_clock;

double ElapsedSeconds(Clock::time_point start) { return std::chrono::duration<double>(Clock::now() - start).count(); }

struct BVHBuildInput {
  LogicalVolume const *logicalVolume{nullptr};
  Vector3D<Precision> *aabbs{nullptr};
  int nDaughters{0};
};

struct BVHBoxProduct {
  LogicalVolume const *logicalVolume{nullptr};
  Vector3D<Precision> *originalAABBs{nullptr};
  int nDaughters{0};
  std::vector<Vector3D<Precision>> boxes;
  std::vector<int> boxToDaughter;
  std::vector<bvhrecipe::Daughter> splitDecisions;
  bool hasSubdivision{false};
  bool fromRecipe{false};
  // Product-local counters are merged after all compute tasks finish. This
  // keeps the split-build phase free of writes to global diagnostics.
  BVHBuildStatistics stats;
};

struct BVHDaughterProduct {
  std::vector<Vector3D<Precision>> boxes;
  int splitAxis{-1};
  int splitM{1};
  bool hasSubdivision{false};
  BVHBuildStatistics stats;
};

struct SubdivisionPrefilterDecision {
  bool trySubdivision{false};
  bool hasFastInflation{false};
  Precision fastInflation{0.};
};

struct BVHDaughterPrefilterProduct {
  SubdivisionPrefilterDecision decision;
  double elapsed{0.};
};

struct RecipeIdentity {
  std::uint64_t fingerprint{0};
};

struct RecipeIdentitySet {
  std::vector<RecipeIdentity> volumes;
  std::uint64_t geometryFingerprint{0};
};

// One independently owned BVH result per logical volume. Workers only write
// their assigned slot; ownership moves into hBVH during serial publication.
struct IndexedBVHConstructionProduct {
  std::unique_ptr<IndexedBVH<BVHBuildReal_t>> bvh;
  double elapsed{0.};
  std::size_t allocatedBytes{0};
};

// One executor job for split-capable modes. The index represents a daughter in
// a mother volume, not the full mother product.
struct BVHDaughterBuildJob {
  std::size_t inputIndex{0};
  int daughterId{0};
};

// One persisted fixed split. Unlike a live daughter job, this carries the
// already selected axis and M and performs no heuristic search.
struct BVHRecipeBuildJob {
  std::size_t inputIndex{0};
  int daughterId{0};
  int axis{0};
  int M{1};
};

// Serial implementation of the same executor contract used by framework and
// std::thread scheduling paths. It keeps the default path behaviorally simple
// and gives Init() one build orchestration flow.
class SerialBVHBuildExecutor final : public BVHBuildExecutor {
public:
  void ParallelFor(std::size_t begin, std::size_t end, std::size_t minGrainSize, BVHBuildRangeTask task,
                   void *context) override
  {
    (void)minGrainSize;
    task(begin, end, context);
  }
};

// Minimal VecGeom-owned scheduler used when no external framework executor is
// installed. It owns worker threads only for the duration of ParallelFor().
class StdThreadBVHBuildExecutor final : public BVHBuildExecutor {
public:
  void ParallelFor(std::size_t begin, std::size_t end, std::size_t minGrainSize, BVHBuildRangeTask task,
                   void *context) override;
};

struct BVHDaughterBuildContext {
  std::vector<BVHBuildInput> const *inputs{nullptr};
  std::vector<BVHDaughterBuildJob> const *jobs{nullptr};
  AABBSubdividerConfig subdivisionConfig;
  // One result slot per job. Workers never append to a shared mother product.
  std::vector<BVHDaughterProduct> *daughterProducts{nullptr};
};

struct BVHDaughterPrefilterContext {
  std::vector<BVHBuildInput> const *inputs{nullptr};
  std::vector<BVHDaughterBuildJob> const *jobs{nullptr};
  AABBSubdividerConfig subdivisionConfig;
  std::vector<std::vector<BVHDaughterPrefilterProduct>> *prefilterProducts{nullptr};
};

struct BVHRecipeBuildContext {
  std::vector<BVHBuildInput> const *inputs{nullptr};
  std::vector<BVHRecipeBuildJob> const *jobs{nullptr};
  std::vector<BVHDaughterProduct> *daughterProducts{nullptr};
};

struct IndexedBVHConstructionContext {
  std::vector<BVHBoxProduct> *boxProducts{nullptr};
  std::vector<std::size_t> const *jobs{nullptr};
  std::vector<IndexedBVHConstructionProduct> *constructionProducts{nullptr};
  AABBSubdividerConfig subdivisionConfig;
};

std::string Lowercase(std::string value)
{
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return value;
}

bool ParseSubdivisionMode(std::string value, AABBSubdivisionMode &mode)
{
  value = Lowercase(value);
  if (value.empty()) return false;
  if (value == "off" || value == "none" || value == "plain" || value == "0") {
    mode = AABBSubdivisionMode::kOff;
    return true;
  }
  if (value == "fast" || value == "1") {
    mode = AABBSubdivisionMode::kFast;
    return true;
  }
  if (value == "balanced" || value == "balance" || value == "2") {
    mode = AABBSubdivisionMode::kBalanced;
    return true;
  }
  if (value == "navigation" || value == "nav" || value == "quality" || value == "aggressive" || value == "3") {
    mode = AABBSubdivisionMode::kNavigation;
    return true;
  }
  VECGEOM_LOG(warning) << "Ignoring unknown " << kSubdivisionModeEnv << " value '" << value << "'";
  return false;
}

void LoadSubdivisionEnvironment()
{
  if (gSubdivisionEnvironmentLoaded) return;
  gSubdivisionEnvironmentLoaded = true;

  if (!gSubdivisionModeExplicit) {
    AABBSubdivisionMode mode;
    if (ParseSubdivisionMode(vecgeom::getenv(kSubdivisionModeEnv), mode)) gSubdivisionConfig.mode = mode;
  }
  if (!gSubdivisionMaxMExplicit) {
    auto const maxM = vecgeom::getenv(kSubdivisionMaxMEnv);
    if (!maxM.empty()) gSubdivisionMaxM = std::atoi(maxM.c_str());
  }
}

void LoadSubdivisionRecipeEnvironment()
{
  if (gSubdivisionRecipeExplicit) return;

  auto const fileName = vecgeom::getenv(kSubdivisionRecipeEnv);
  if (fileName == gSubdivisionRecipeEnvironmentFile) return;
  gSubdivisionRecipeCreateFile.clear();
  gLoadedRecipe.reset();
  if (fileName.empty()) {
    gSubdivisionRecipeEnvironmentFile.clear();
    return;
  }

  std::error_code statusError;
  bool const exists = std::filesystem::exists(fileName, statusError);
  if (statusError) {
    // Keep failed paths retryable so a transient filesystem error does not
    // disable the environment recipe for the rest of the process.
    gSubdivisionRecipeEnvironmentFile.clear();
    VECGEOM_LOG(warning) << "Cannot inspect " << kSubdivisionRecipeEnv << "='" << fileName
                         << "': " << statusError.message();
    return;
  }
  if (!exists) {
    // Defer creation until every BVH has been built and published successfully.
    gSubdivisionRecipeEnvironmentFile = fileName;
    gSubdivisionRecipeCreateFile      = fileName;
    return;
  }

  auto recipe = std::make_unique<bvhrecipe::Recipe>();
  std::string error;
  if (!bvhrecipe::Read(fileName, *recipe, error)) {
    // A repaired or atomically replaced file at the same path must be observed
    // by the next Init(), so only successful reads are cached.
    gSubdivisionRecipeEnvironmentFile.clear();
    VECGEOM_LOG(warning) << "Ignoring " << kSubdivisionRecipeEnv << "='" << fileName << "': " << error;
    return;
  }
  gSubdivisionRecipeEnvironmentFile = fileName;
  gLoadedRecipe                     = std::move(recipe);
}

AABBSubdividerConfig CurrentSubdivisionConfig()
{
  LoadSubdivisionEnvironment();
  AABBSubdividerConfig cfg = gSubdivisionConfig;
  if (gSubdivisionMaxM > 0) {
    cfg.fastMaxM     = std::min(gSubdivisionMaxM, kMaximumAABBSubdivisionM);
    cfg.balancedMaxM = cfg.fastMaxM;
    cfg.autoMaxM     = cfg.fastMaxM;
  }
  cfg.fastMaxM     = std::clamp(cfg.fastMaxM, 1, kMaximumAABBSubdivisionM);
  cfg.balancedMaxM = std::clamp(cfg.balancedMaxM, 1, kMaximumAABBSubdivisionM);
  cfg.autoMaxM     = std::clamp(cfg.autoMaxM, 1, kMaximumAABBSubdivisionM);
  return cfg;
}

// Stable FNV-1a hashing with explicit byte order. Persisted identities must not
// depend on host endianness or pointer values.
class RecipeHash {
private:
  std::uint64_t fValue = UINT64_C(14695981039346656037);

  void AddByte(unsigned char value)
  {
    fValue ^= value;
    fValue *= UINT64_C(1099511628211);
  }

public:
  void AddUnsigned(std::uint64_t value, std::size_t bytes = sizeof(value))
  {
    for (std::size_t i = 0; i < bytes; ++i) {
      AddByte(static_cast<unsigned char>(value & UINT64_C(0xff)));
      value >>= 8;
    }
  }

  void AddString(std::string const &value)
  {
    AddUnsigned(value.size());
    for (unsigned char c : value)
      AddByte(c);
  }

  void AddPrecision(Precision value)
  {
    static_assert(sizeof(Precision) <= sizeof(std::uint64_t), "unsupported Precision size");
    std::uint64_t bits = 0;
    std::memcpy(&bits, &value, sizeof(value));
    AddUnsigned(bits, sizeof(value));
  }

  std::uint64_t Value() const { return fValue; }
};

std::uint64_t ShapeFingerprint(VUnplacedVolume const *unplaced)
{
  // The file stores axis/M, not box coordinates. All slice boxes are calculated
  // again from the current solid, even when a saved split choice is reused.
  RecipeHash hash;
  hash.AddUnsigned(static_cast<unsigned char>(unplaced->GetType()), sizeof(unsigned char));
  hash.AddUnsigned(static_cast<unsigned>(unplaced->MemorySize()), sizeof(unsigned));
  return hash.Value();
}

std::uint64_t DaughterFingerprint(VPlacedVolume const *daughter, Vector3D<Precision> const &aabbLower,
                                  Vector3D<Precision> const &aabbUpper, std::uint64_t shapeFingerprint)
{
  RecipeHash hash;
  hash.AddString(daughter->GetLabel());
  hash.AddString(daughter->GetLogicalVolume()->GetLabel());
  hash.AddUnsigned(shapeFingerprint);

  auto const *transform = daughter->GetTransformation();
  for (int axis = 0; axis < 3; ++axis)
    hash.AddPrecision(transform->Translation(axis));
  for (int element = 0; element < 9; ++element)
    hash.AddPrecision(transform->Rotation(element));
  for (int axis = 0; axis < 3; ++axis) {
    hash.AddPrecision(aabbLower[axis]);
    hash.AddPrecision(aabbUpper[axis]);
  }
  return hash.Value();
}

RecipeIdentitySet ComputeRecipeIdentities(std::vector<BVHBuildInput> const &inputs)
{
  RecipeIdentitySet identities;
  identities.volumes.resize(inputs.size());
  std::unordered_map<VUnplacedVolume const *, std::uint64_t> shapeFingerprints;
  std::vector<std::pair<std::string, std::uint64_t>> geometryEntries;

  auto getShapeFingerprint = [&shapeFingerprints](VUnplacedVolume const *unplaced) {
    auto inserted = shapeFingerprints.emplace(unplaced, 0);
    if (inserted.second) inserted.first->second = ShapeFingerprint(unplaced);
    return inserted.first->second;
  };

  for (std::size_t inputIndex = 0; inputIndex < inputs.size(); ++inputIndex) {
    auto const &input = inputs[inputIndex];
    auto &identity    = identities.volumes[inputIndex];
    RecipeHash volumeHash;
    volumeHash.AddString(input.logicalVolume->GetLabel());
    volumeHash.AddUnsigned(getShapeFingerprint(input.logicalVolume->GetUnplacedVolume()));
    volumeHash.AddUnsigned(static_cast<unsigned>(input.nDaughters), sizeof(unsigned));
    for (int daughterId = 0; daughterId < input.nDaughters; ++daughterId) {
      auto const *daughter = input.logicalVolume->GetDaughters()[daughterId];
      auto const fingerprint =
          DaughterFingerprint(daughter, input.aabbs[2 * daughterId], input.aabbs[2 * daughterId + 1],
                              getShapeFingerprint(daughter->GetUnplacedVolume()));
      volumeHash.AddUnsigned(fingerprint);
    }
    identity.fingerprint = volumeHash.Value();
    if (input.nDaughters > 0) geometryEntries.emplace_back(input.logicalVolume->GetLabel(), identity.fingerprint);
  }

  std::sort(geometryEntries.begin(), geometryEntries.end());
  RecipeHash geometryHash;
  for (auto const &entry : geometryEntries) {
    geometryHash.AddString(entry.first);
    geometryHash.AddUnsigned(entry.second);
  }
  identities.geometryFingerprint = geometryHash.Value();
  return identities;
}

Precision SurfaceArea(Vector3D<Precision> const &lower, Vector3D<Precision> const &upper)
{
  auto extent = upper - lower;
  if (extent[0] <= 0 || extent[1] <= 0 || extent[2] <= 0) return 0;
  return 2 * (extent[0] * extent[1] + extent[1] * extent[2] + extent[2] * extent[0]);
}

SubdivisionPrefilterDecision EvaluateSubdivisionPrefilter(VPlacedVolume const *daughter,
                                                          Vector3D<Precision> const &aabbLower,
                                                          Vector3D<Precision> const &aabbUpper,
                                                          AABBSubdividerConfig const &cfg)
{
  SubdivisionPrefilterDecision decision;
  if (cfg.mode == AABBSubdivisionMode::kOff) return decision;
  if (!daughter || daughter->IsAssembly()) return decision;
  const Precision globalAABBSA = SurfaceArea(aabbLower, aabbUpper);
  if (globalAABBSA <= 0.) return decision;

  if (cfg.mode == AABBSubdivisionMode::kNavigation) {
    decision.trySubdivision = true;
    return decision;
  }

  Vector3D<Precision> extentLower, extentUpper;
  daughter->GetUnplacedVolume()->Extent(extentLower, extentUpper);

  const Precision localExtentSA = SurfaceArea(extentLower, extentUpper);
  if (localExtentSA <= 0.) return decision;
  decision.hasFastInflation = true;
  decision.fastInflation    = globalAABBSA / localExtentSA;
  const Precision minInflation =
      cfg.mode == AABBSubdivisionMode::kBalanced ? cfg.balancedMinAABBInflation : cfg.fastMinAABBInflation;
  decision.trySubdivision = decision.fastInflation > minInflation;
  return decision;
}

bool AcceptSubdivision(AABBSubdivisionResult const &result)
{
  return result.chosenM > 1 && !result.slabs.empty() && result.globalSA > 0 && result.finalSA < result.globalSA;
}

BVHBuildExecutor &CurrentBuildExecutor()
{
  // External frameworks get priority to avoid oversubscribing applications that
  // already control geometry-initialization parallelism.
  if (gBuildExecutor) return *gBuildExecutor;
  static StdThreadBVHBuildExecutor stdThreadExecutor;
  static SerialBVHBuildExecutor serialExecutor;
  return gBuildThreadCount == 1 ? static_cast<BVHBuildExecutor &>(serialExecutor)
                                : static_cast<BVHBuildExecutor &>(stdThreadExecutor);
}

std::size_t ResolveBuildThreadCount(std::size_t jobs, std::size_t minGrainSize)
{
  if (jobs == 0) return 1;

  unsigned requested = gBuildThreadCount;
  if (requested == BVHManager::kCPUCount) {
    requested = std::thread::hardware_concurrency();
    if (requested == 0) requested = 1;
  }

  const std::size_t grain      = std::max<std::size_t>(1, minGrainSize);
  const std::size_t usefulJobs = (jobs + grain - 1) / grain;
  // Never spawn more workers than useful chunks. This keeps small geometries
  // and explicit large caps from paying avoidable thread overhead.
  return std::max<std::size_t>(1, std::min<std::size_t>(requested, usefulJobs));
}

void StdThreadBVHBuildExecutor::ParallelFor(std::size_t begin, std::size_t end, std::size_t minGrainSize,
                                            BVHBuildRangeTask task, void *context)
{
  if (begin >= end) return;

  const std::size_t threadCount = ResolveBuildThreadCount(end - begin, minGrainSize);
  if (threadCount <= 1) {
    task(begin, end, context);
    return;
  }

  // Workers claim fixed-size ranges atomically. The VecGeom task owns the
  // meaning of range indices and writes to disjoint product slots.
  const std::size_t grain = std::max<std::size_t>(1, minGrainSize);
  std::atomic<std::size_t> next{begin};
  std::atomic<bool> failed{false};
  std::exception_ptr firstException;
  std::mutex exceptionMutex;

  auto worker = [&]() {
    while (!failed.load(std::memory_order_relaxed)) {
      const std::size_t rangeBegin = next.fetch_add(grain, std::memory_order_relaxed);
      if (rangeBegin >= end) break;
      const std::size_t rangeEnd = std::min(end, rangeBegin + grain);
      try {
        task(rangeBegin, rangeEnd, context);
      } catch (...) {
        if (!failed.exchange(true, std::memory_order_relaxed)) {
          std::lock_guard<std::mutex> lock(exceptionMutex);
          firstException = std::current_exception();
        }
        break;
      }
    }
  };

  std::vector<std::thread> threads;
  threads.reserve(threadCount);
  try {
    for (std::size_t i = 0; i < threadCount; ++i) {
      threads.emplace_back(worker);
    }
  } catch (...) {
    // Destroying a joinable std::thread terminates the process. Stop range
    // assignment and join workers that started before creation failed.
    failed.store(true, std::memory_order_relaxed);
    for (auto &thread : threads)
      thread.join();
    throw;
  }
  for (auto &thread : threads) {
    thread.join();
  }

  // Re-throw the first worker failure after all threads have joined so Init()
  // has the same synchronous failure semantics for every executor.
  if (firstException) std::rethrow_exception(firstException);
}

void AccumulateProductStats(BVHBuildStatistics &stats, BVHBuildStatistics const &productStats)
{
  stats.candidateSelectionElapsed += productStats.candidateSelectionElapsed;
  stats.subdivisionElapsed += productStats.subdivisionElapsed;
  stats.splitCandidates += productStats.splitCandidates;
  stats.balancedSingleAxisCandidates += productStats.balancedSingleAxisCandidates;
  stats.balancedMultiAxisCandidates += productStats.balancedMultiAxisCandidates;
  stats.acceptedSplits += productStats.acceptedSplits;
  stats.splitBoxes += productStats.splitBoxes;
  stats.extraBoxes += productStats.extraBoxes;
  stats.costEvaluatedLogicalVolumes += productStats.costEvaluatedLogicalVolumes;
  stats.costRejectedLogicalVolumes += productStats.costRejectedLogicalVolumes;
  stats.distanceCostRejectedLogicalVolumes += productStats.distanceCostRejectedLogicalVolumes;
  stats.costRejectedSplits += productStats.costRejectedSplits;
  stats.costRejectedExtraBoxes += productStats.costRejectedExtraBoxes;
  stats.fastMotherSkippedDaughters += productStats.fastMotherSkippedDaughters;
  stats.fastInflationSamples += productStats.fastInflationSamples;
  stats.fastInflationSum += productStats.fastInflationSum;
  stats.fastInflationMax = std::max(stats.fastInflationMax, productStats.fastInflationMax);
  stats.fastInflationAbove1p5 += productStats.fastInflationAbove1p5;
  stats.fastInflationAbove2 += productStats.fastInflationAbove2;
  stats.fastInflationAbove3 += productStats.fastInflationAbove3;
  stats.fastInflationAbove5 += productStats.fastInflationAbove5;
  stats.fastInflationAbove8 += productStats.fastInflationAbove8;
  stats.fastInflationAbove13 += productStats.fastInflationAbove13;
}

bool FastMotherWithinSizeLimit(BVHBuildInput const &input, AABBSubdividerConfig const &config)
{
  return config.mode != AABBSubdivisionMode::kFast || config.fastMaxDaughtersPerMother <= 0 ||
         input.nDaughters <= config.fastMaxDaughtersPerMother;
}

bool FastMotherHasEnoughCandidates(std::size_t candidates, int nDaughters, AABBSubdividerConfig const &config)
{
  return config.mode != AABBSubdivisionMode::kFast || config.fastMinCandidateFractionPerMother <= 0. ||
         Precision(candidates) >= config.fastMinCandidateFractionPerMother * Precision(nDaughters);
}

void RecordFastInflation(BVHBuildStatistics &stats, Precision ratio)
{
  ++stats.fastInflationSamples;
  stats.fastInflationSum += ratio;
  stats.fastInflationMax = std::max(stats.fastInflationMax, double(ratio));
  if (ratio > Precision(1.5)) ++stats.fastInflationAbove1p5;
  if (ratio > Precision(2.0)) ++stats.fastInflationAbove2;
  if (ratio > Precision(3.0)) ++stats.fastInflationAbove3;
  if (ratio > Precision(5.0)) ++stats.fastInflationAbove5;
  if (ratio > Precision(8.0)) ++stats.fastInflationAbove8;
  if (ratio > Precision(13.0)) ++stats.fastInflationAbove13;
}

void InitializeProduct(BVHBoxProduct &product, BVHBuildInput const &input)
{
  product.logicalVolume = input.logicalVolume;
  product.originalAABBs = input.aabbs;
  product.nDaughters    = input.nDaughters;
}

void AppendOriginalDaughterBox(BVHBoxProduct &product, int daughterId);
void AppendDaughterProduct(BVHBoxProduct &product, int daughterId, BVHDaughterProduct const &daughterProduct);

bool RecipeMetadataCompatible(bvhrecipe::Recipe const &recipe, std::string &reason)
{
  if (recipe.precisionBytes != static_cast<int>(sizeof(Precision))) {
    reason = "Precision size differs";
    return false;
  }
  if (recipe.bvhRealBytes != static_cast<int>(sizeof(BVHBuildReal_t))) {
    reason = "BVH real size differs";
    return false;
  }
  return true;
}

bool RecipeVolumeMatches(bvhrecipe::Volume const &recipeVolume, BVHBuildInput const &input,
                         RecipeIdentity const &identity)
{
  return recipeVolume.label == input.logicalVolume->GetLabel() && recipeVolume.fingerprint == identity.fingerprint &&
         recipeVolume.daughterCount == static_cast<std::size_t>(input.nDaughters);
}

BVHDaughterProduct ReconstructRecipeSubdivision(BVHBuildInput const &input, BVHRecipeBuildJob const &job)
{
  BVHDaughterProduct product;
  auto const *daughter = input.logicalVolume->GetDaughters()[job.daughterId];
  SlabBox baseline;
  baseline.min                     = input.aabbs[2 * job.daughterId];
  baseline.max                     = input.aabbs[2 * job.daughterId + 1];
  auto const start                 = Clock::now();
  auto slabs                       = AABBSubdivider::Instance().SubdivideAlongAxis(daughter, job.axis, job.M, baseline);
  product.stats.subdivisionElapsed = ElapsedSeconds(start);
  if (slabs.empty()) return product;

  product.boxes.reserve(2 * slabs.size());
  for (auto const &slab : slabs) {
    product.boxes.push_back(slab.min);
    product.boxes.push_back(slab.max);
  }
  product.splitAxis      = job.axis;
  product.splitM         = job.M;
  product.hasSubdivision = true;
  ++product.stats.acceptedSplits;
  product.stats.splitBoxes += slabs.size();
  product.stats.extraBoxes += slabs.size() - 1;
  return product;
}

void ReconstructRecipeSubdivisionRange(std::size_t begin, std::size_t end, void *context)
{
  auto &buildContext = *static_cast<BVHRecipeBuildContext *>(context);
  for (std::size_t index = begin; index < end; ++index) {
    auto const &job                         = (*buildContext.jobs)[index];
    (*buildContext.daughterProducts)[index] = ReconstructRecipeSubdivision((*buildContext.inputs)[job.inputIndex], job);
  }
}

bool AssembleRecipeVolume(BVHBoxProduct &product, BVHBuildInput const &input,
                          std::vector<BVHRecipeBuildJob> const &jobs,
                          std::vector<BVHDaughterProduct> const &daughterProducts, std::size_t first, std::size_t last,
                          BVHBuildStatistics &buildStats)
{
  for (std::size_t jobIndex = first; jobIndex < last; ++jobIndex) {
    if (!daughterProducts[jobIndex].hasSubdivision) return false;
  }
  if (first == last) return true;

  std::size_t primitiveCount = input.nDaughters;
  for (std::size_t jobIndex = first; jobIndex < last; ++jobIndex)
    primitiveCount += daughterProducts[jobIndex].stats.extraBoxes;
  product.boxes.reserve(2 * primitiveCount);
  product.boxToDaughter.reserve(primitiveCount);
  product.splitDecisions.reserve(last - first);

  std::size_t jobIndex = first;
  for (int daughterId = 0; daughterId < input.nDaughters; ++daughterId) {
    if (jobIndex < last && jobs[jobIndex].daughterId == daughterId) {
      auto const boxCount = daughterProducts[jobIndex].boxes.size() / 2;
      ++buildStats.recipeSplitDecisionsUsed;
      buildStats.recipeSplitBoxesBuilt += boxCount;
      AppendDaughterProduct(product, daughterId, daughterProducts[jobIndex]);
      ++jobIndex;
    } else {
      AppendOriginalDaughterBox(product, daughterId);
    }
  }
  return true;
}

std::vector<bool> ApplyLoadedRecipe(std::vector<BVHBuildInput> const &inputs, std::vector<BVHBoxProduct> &products,
                                    BVHBuildStatistics &buildStats)
{
  std::vector<bool> liveBuild(inputs.size(), true);
  if (!gLoadedRecipe) return liveBuild;

  // A recipe is a one-build input. Environment selection reloads it if a later
  // Init() explicitly rebuilds the geometry.
  auto recipe = std::move(gLoadedRecipe);
  if (!gSubdivisionRecipeExplicit) gSubdivisionRecipeEnvironmentFile.clear();

  auto const validationStart     = Clock::now();
  buildStats.recipeEntriesLoaded = recipe->volumes.size();
  auto countLiveMothers          = [&]() {
    return std::count_if(inputs.begin(), inputs.end(), [](BVHBuildInput const &input) { return input.nDaughters > 0; });
  };

  std::string incompatibility;
  if (!RecipeMetadataCompatible(*recipe, incompatibility)) {
    buildStats.recipeEntriesRejected   = recipe->volumes.size();
    buildStats.recipeEntriesLiveBuilt  = countLiveMothers();
    buildStats.recipeValidationElapsed = ElapsedSeconds(validationStart);
    VECGEOM_LOG(warning) << "Ignoring BVH AABB subdivision recipe: " << incompatibility;
    return liveBuild;
  }
  if (recipe->vecgeomVersion != vecgeom_version) {
    VECGEOM_LOG(info) << "Reusing BVH AABB split decisions produced by VecGeom " << recipe->vecgeomVersion << " with "
                      << vecgeom_version << "; slab boxes are rebuilt by the current implementation";
  }

  auto const identities        = ComputeRecipeIdentities(inputs);
  auto const *world            = GeoManager::Instance().GetWorld();
  std::string const worldLabel = world ? world->GetLabel() : std::string{};
  if (recipe->worldLabel != worldLabel || recipe->geometryFingerprint != identities.geometryFingerprint) {
    VECGEOM_LOG(warning) << "BVH AABB subdivision recipe geometry identity differs; validating entries individually";
  }

  std::vector<bool> used(recipe->volumes.size(), false);
  std::vector<std::size_t> matchedRecipe(inputs.size(), recipe->volumes.size());
  std::unordered_map<int, std::size_t> recipeByLogicalId;
  std::unordered_multimap<std::uint64_t, std::size_t> recipeByFingerprint;
  for (std::size_t recipeIndex = 0; recipeIndex < recipe->volumes.size(); ++recipeIndex) {
    auto const &volume = recipe->volumes[recipeIndex];
    recipeByLogicalId.emplace(volume.logicalId, recipeIndex);
    recipeByFingerprint.emplace(volume.fingerprint, recipeIndex);
  }

  for (std::size_t inputIndex = 0; inputIndex < inputs.size(); ++inputIndex) {
    auto const &input = inputs[inputIndex];
    if (input.nDaughters == 0) {
      liveBuild[inputIndex] = false;
      continue;
    }

    std::size_t match  = recipe->volumes.size();
    auto const idMatch = recipeByLogicalId.find(static_cast<int>(input.logicalVolume->id()));
    if (idMatch != recipeByLogicalId.end() && !used[idMatch->second] &&
        RecipeVolumeMatches(recipe->volumes[idMatch->second], input, identities.volumes[inputIndex])) {
      match = idMatch->second;
    }
    if (match == recipe->volumes.size()) {
      auto const candidates = recipeByFingerprint.equal_range(identities.volumes[inputIndex].fingerprint);
      for (auto candidate = candidates.first; candidate != candidates.second; ++candidate) {
        if (!used[candidate->second] &&
            RecipeVolumeMatches(recipe->volumes[candidate->second], input, identities.volumes[inputIndex])) {
          match = candidate->second;
          break;
        }
      }
    }

    if (match == recipe->volumes.size()) {
      ++buildStats.recipeEntriesMissing;
      continue;
    }

    used[match]               = true;
    matchedRecipe[inputIndex] = match;
    liveBuild[inputIndex]     = false;
  }

  const double matchingElapsed = ElapsedSeconds(validationStart);

  // Recipe records are canonicalized by daughter id, so jobs are already
  // grouped by mother and can be assembled without an additional sort.
  std::vector<std::size_t> firstRecipeJob(inputs.size() + 1, 0);
  std::vector<BVHRecipeBuildJob> recipeJobs;
  for (std::size_t inputIndex = 0; inputIndex < inputs.size(); ++inputIndex) {
    firstRecipeJob[inputIndex] = recipeJobs.size();
    auto const recipeIndex     = matchedRecipe[inputIndex];
    if (recipeIndex == recipe->volumes.size()) continue;
    for (auto const &decision : recipe->volumes[recipeIndex].splitDaughters)
      recipeJobs.push_back({inputIndex, decision.daughterId, decision.axis, decision.M});
  }
  firstRecipeJob[inputs.size()] = recipeJobs.size();

  std::vector<BVHDaughterProduct> recipeProducts(recipeJobs.size());
  BVHRecipeBuildContext recipeContext{&inputs, &recipeJobs, &recipeProducts};
  auto const reconstructionStart = Clock::now();
  CurrentBuildExecutor().ParallelFor(0, recipeJobs.size(), 1, &ReconstructRecipeSubdivisionRange, &recipeContext);
  buildStats.recipeRebuildWallElapsed = ElapsedSeconds(reconstructionStart);
  for (auto &product : recipeProducts) {
    buildStats.recipeRebuildElapsed += product.stats.subdivisionElapsed;
    product.stats.subdivisionElapsed = 0.;
  }

  auto const assemblyStart = Clock::now();
  for (std::size_t inputIndex = 0; inputIndex < inputs.size(); ++inputIndex) {
    auto const recipeIndex = matchedRecipe[inputIndex];
    if (recipeIndex == recipe->volumes.size()) continue;

    bool const reconstructed =
        AssembleRecipeVolume(products[inputIndex], inputs[inputIndex], recipeJobs, recipeProducts,
                             firstRecipeJob[inputIndex], firstRecipeJob[inputIndex + 1], buildStats);
    if (reconstructed) {
      products[inputIndex].fromRecipe = true;
      ++buildStats.recipeEntriesUsed;
      continue;
    }

    used[recipeIndex]                   = false;
    liveBuild[inputIndex]               = true;
    matchedRecipe[inputIndex]           = recipe->volumes.size();
    products[inputIndex].boxes          = {};
    products[inputIndex].boxToDaughter  = {};
    products[inputIndex].splitDecisions = {};
    products[inputIndex].hasSubdivision = false;
    products[inputIndex].fromRecipe     = false;
    products[inputIndex].stats          = {};
    VECGEOM_LOG(warning) << "Cannot reconstruct persisted BVH AABB subdivision for logical volume '"
                         << inputs[inputIndex].logicalVolume->GetLabel() << "'; using the live fallback mode";
  }

  buildStats.recipeEntriesRejected = std::count(used.begin(), used.end(), false);
  for (std::size_t inputIndex = 0; inputIndex < inputs.size(); ++inputIndex) {
    if (liveBuild[inputIndex] && inputs[inputIndex].nDaughters > 0) ++buildStats.recipeEntriesLiveBuilt;
  }
  buildStats.recipeValidationElapsed = matchingElapsed + ElapsedSeconds(assemblyStart);
  return liveBuild;
}

std::vector<BVHBuildInput> SnapshotBuildInputs(std::vector<LogicalVolume const *> const &logicalVolumes,
                                               BVHBuildStatistics &stats)
{
  std::vector<BVHBuildInput> inputs;
  inputs.reserve(logicalVolumes.size());
  for (auto logicalVolume : logicalVolumes) {
    ++stats.logicalVolumes;

    int nDaughters{0};
    auto *aabbs = ABBoxManager<Precision>::Instance().GetABBoxes(logicalVolume, nDaughters);
    if (nDaughters > 0) {
      ++stats.logicalVolumesWithBoxes;
      stats.originalBoxes += static_cast<std::size_t>(nDaughters);
    }
    inputs.push_back({logicalVolume, aabbs, nDaughters});
  }
  return inputs;
}

BVHDaughterProduct ComputeSubdivisionProduct(BVHBuildInput const &input, int daughterId,
                                             AABBSubdividerConfig const &subdivisionConfig)
{
  BVHDaughterProduct product;
  auto &subdivider = AABBSubdivider::Instance();

  VPlacedVolume const *daughter = input.logicalVolume->GetDaughters()[daughterId];
  ++product.stats.splitCandidates;
  auto subdivisionStart = Clock::now();
  auto result           = subdivider.SubdivideAuto(daughter, subdivisionConfig);
  product.stats.subdivisionElapsed += ElapsedSeconds(subdivisionStart);
  if (subdivisionConfig.mode == AABBSubdivisionMode::kBalanced) {
    if (result.usedSingleAxisSearch)
      ++product.stats.balancedSingleAxisCandidates;
    else
      ++product.stats.balancedMultiAxisCandidates;
  }
  if (AcceptSubdivision(result)) {
    product.boxes.reserve(2 * result.slabs.size());
    for (auto const &slab : result.slabs) {
      product.boxes.push_back(slab.min);
      product.boxes.push_back(slab.max);
    }
    product.splitAxis      = result.chosenAxis;
    product.splitM         = result.chosenM;
    product.hasSubdivision = true;
    ++product.stats.acceptedSplits;
    product.stats.splitBoxes += result.slabs.size();
    product.stats.extraBoxes += result.slabs.size() > 0 ? result.slabs.size() - 1 : 0;
  }
  return product;
}

void AppendOriginalDaughterBox(BVHBoxProduct &product, int daughterId)
{
  product.boxes.push_back(product.originalAABBs[2 * daughterId]);
  product.boxes.push_back(product.originalAABBs[2 * daughterId + 1]);
  product.boxToDaughter.push_back(daughterId);
}

void AppendDaughterProduct(BVHBoxProduct &product, int daughterId, BVHDaughterProduct const &daughterProduct)
{
  AccumulateProductStats(product.stats, daughterProduct.stats);
  product.hasSubdivision = product.hasSubdivision || daughterProduct.hasSubdivision;

  if (!daughterProduct.hasSubdivision) {
    AppendOriginalDaughterBox(product, daughterId);
    return;
  }
  const auto nBoxes = daughterProduct.boxes.size() / 2;
  product.boxes.insert(product.boxes.end(), daughterProduct.boxes.begin(), daughterProduct.boxes.end());
  for (std::size_t i = 0; i < nBoxes; ++i)
    product.boxToDaughter.push_back(daughterId);

  bvhrecipe::Daughter decision;
  decision.daughterId = daughterId;
  decision.axis       = daughterProduct.splitAxis;
  decision.M          = daughterProduct.splitM;
  product.splitDecisions.push_back(decision);
}

void ComputeDaughterProductRange(std::size_t begin, std::size_t end, void *context)
{
  auto &buildContext = *static_cast<BVHDaughterBuildContext *>(context);
  for (std::size_t index = begin; index < end; ++index) {
    auto const &job                         = (*buildContext.jobs)[index];
    (*buildContext.daughterProducts)[index] = ComputeSubdivisionProduct((*buildContext.inputs)[job.inputIndex],
                                                                        job.daughterId, buildContext.subdivisionConfig);
  }
}

void ComputeDaughterPrefilterRange(std::size_t begin, std::size_t end, void *context)
{
  auto &buildContext = *static_cast<BVHDaughterPrefilterContext *>(context);
  for (std::size_t index = begin; index < end; ++index) {
    auto const &job      = (*buildContext.jobs)[index];
    auto const &input    = (*buildContext.inputs)[job.inputIndex];
    auto const *daughter = input.logicalVolume->GetDaughters()[job.daughterId];
    auto &product        = (*buildContext.prefilterProducts)[job.inputIndex][job.daughterId];
    auto start           = Clock::now();
    product.decision     = EvaluateSubdivisionPrefilter(
        daughter, input.aabbs[2 * job.daughterId], input.aabbs[2 * job.daughterId + 1], buildContext.subdivisionConfig);
    product.elapsed = ElapsedSeconds(start);
  }
}

std::vector<BVHBoxProduct> BuildBoxProducts(std::vector<BVHBuildInput> const &inputs,
                                            AABBSubdividerConfig const &subdivisionConfig,
                                            BVHBuildStatistics &buildStats)
{
  // Build in two deterministic phases. The first phase evaluates only the
  // cheap daughter prefilter and applies the fast-mode policy per mother. The
  // second creates clipped-extent jobs only for admitted candidates. Mothers
  // with no accepted split remain identity BVHs and never allocate mapped box
  // products proportional to their full daughter count.
  std::vector<BVHBoxProduct> products(inputs.size());
  for (std::size_t index = 0; index < inputs.size(); ++index)
    InitializeProduct(products[index], inputs[index]);
  auto const liveBuild = ApplyLoadedRecipe(inputs, products, buildStats);
  if (subdivisionConfig.mode == AABBSubdivisionMode::kOff) return products;

  std::vector<std::vector<BVHDaughterPrefilterProduct>> prefilterProducts(inputs.size());
  std::vector<BVHDaughterBuildJob> prefilterJobs;

  for (std::size_t inputIndex = 0; inputIndex < inputs.size(); ++inputIndex) {
    if (!liveBuild[inputIndex]) continue;
    if (!FastMotherWithinSizeLimit(inputs[inputIndex], subdivisionConfig)) {
      products[inputIndex].stats.fastMotherSkippedDaughters += inputs[inputIndex].nDaughters;
      continue;
    }
    prefilterProducts[inputIndex].resize(inputs[inputIndex].nDaughters);
    for (int daughterId = 0; daughterId < inputs[inputIndex].nDaughters; ++daughterId)
      prefilterJobs.push_back({inputIndex, daughterId});
  }

  BVHDaughterPrefilterContext prefilterContext{&inputs, &prefilterJobs, subdivisionConfig, &prefilterProducts};
  auto candidateSelectionStart = Clock::now();
  CurrentBuildExecutor().ParallelFor(0, prefilterJobs.size(), 1, &ComputeDaughterPrefilterRange, &prefilterContext);
  buildStats.candidateSelectionWallElapsed = ElapsedSeconds(candidateSelectionStart);

  auto candidatePreparationStart = Clock::now();
  std::vector<std::size_t> firstCandidateJob(inputs.size() + 1, 0);
  std::vector<BVHDaughterBuildJob> candidateJobs;
  for (std::size_t inputIndex = 0; inputIndex < inputs.size(); ++inputIndex) {
    firstCandidateJob[inputIndex] = candidateJobs.size();
    auto const &inputPrefilters   = prefilterProducts[inputIndex];
    if (inputPrefilters.empty()) continue;

    std::size_t candidates = 0;
    for (auto const &prefilter : inputPrefilters) {
      products[inputIndex].stats.candidateSelectionElapsed += prefilter.elapsed;
      if (prefilter.decision.hasFastInflation)
        RecordFastInflation(products[inputIndex].stats, prefilter.decision.fastInflation);
      if (prefilter.decision.trySubdivision) ++candidates;
    }
    if (!FastMotherHasEnoughCandidates(candidates, inputs[inputIndex].nDaughters, subdivisionConfig)) {
      products[inputIndex].stats.fastMotherSkippedDaughters += inputs[inputIndex].nDaughters;
      continue;
    }

    for (int daughterId = 0; daughterId < inputs[inputIndex].nDaughters; ++daughterId) {
      if (inputPrefilters[daughterId].decision.trySubdivision) candidateJobs.push_back({inputIndex, daughterId});
    }
  }
  firstCandidateJob[inputs.size()]       = candidateJobs.size();
  buildStats.candidatePreparationElapsed = ElapsedSeconds(candidatePreparationStart);

  // Each worker writes one daughter result slot. The final mother-product arrays
  // are assembled below in stable daughter order, preserving deterministic BVHs.
  std::vector<BVHDaughterProduct> daughterProducts(candidateJobs.size());
  BVHDaughterBuildContext buildContext{&inputs, &candidateJobs, subdivisionConfig, &daughterProducts};
  auto subdivisionStart = Clock::now();
  CurrentBuildExecutor().ParallelFor(0, candidateJobs.size(), 1, &ComputeDaughterProductRange, &buildContext);
  buildStats.subdivisionWallElapsed = ElapsedSeconds(subdivisionStart);

  auto boxAssemblyStart = Clock::now();
  for (std::size_t inputIndex = 0; inputIndex < inputs.size(); ++inputIndex) {
    if (!liveBuild[inputIndex]) continue;
    const std::size_t first   = firstCandidateJob[inputIndex];
    const std::size_t last    = firstCandidateJob[inputIndex + 1];
    const bool hasSubdivision = std::any_of(daughterProducts.begin() + first, daughterProducts.begin() + last,
                                            [](BVHDaughterProduct const &product) { return product.hasSubdivision; });
    if (!hasSubdivision) {
      for (std::size_t jobIndex = first; jobIndex < last; ++jobIndex)
        AccumulateProductStats(products[inputIndex].stats, daughterProducts[jobIndex].stats);
      continue;
    }

    std::size_t primitiveCount = inputs[inputIndex].nDaughters;
    for (std::size_t jobIndex = first; jobIndex < last; ++jobIndex)
      primitiveCount += daughterProducts[jobIndex].stats.extraBoxes;
    products[inputIndex].boxes.reserve(2 * primitiveCount);
    products[inputIndex].boxToDaughter.reserve(primitiveCount);
    products[inputIndex].splitDecisions.reserve(
        std::count_if(daughterProducts.begin() + first, daughterProducts.begin() + last,
                      [](BVHDaughterProduct const &product) { return product.hasSubdivision; }));

    std::size_t jobIndex = first;
    for (int daughterId = 0; daughterId < inputs[inputIndex].nDaughters; ++daughterId) {
      if (jobIndex < last && candidateJobs[jobIndex].daughterId == daughterId) {
        AppendDaughterProduct(products[inputIndex], daughterId, daughterProducts[jobIndex]);
        ++jobIndex;
      } else {
        AppendOriginalDaughterBox(products[inputIndex], daughterId);
      }
    }
  }
  buildStats.boxAssemblyElapsed = ElapsedSeconds(boxAssemblyStart);
  return products;
}

Precision MinimumBVHCostReduction(AABBSubdividerConfig const &config)
{
  if (config.mode == AABBSubdivisionMode::kNavigation) return config.navigationMinCostReduction;
  if (config.mode == AABBSubdivisionMode::kBalanced) return config.balancedMinCostReduction;
  return config.fastMinCostReduction;
}

double EstimatedBVHCost(IndexedBVH<BVHBuildReal_t> const &bvh, AABBSubdividerConfig const &config)
{
  double nodes = 0., boxes = 0., candidates = 0.;
  bvh.EstimateTraversalWork(nodes, boxes, candidates);
  return static_cast<double>(config.nodeCost) * nodes + static_cast<double>(config.boxCost) * boxes +
         static_cast<double>(config.solidCost) * candidates;
}

struct SplitBVHCostEvaluation {
  std::unique_ptr<IndexedBVH<BVHBuildReal_t>> originalBVH;
  bool distanceAccepted{true};
};

SplitBVHCostEvaluation EvaluateSplitBVHCost(IndexedBVH<BVHBuildReal_t> const &splitBVH,
                                            Vector3D<Precision> *originalAABBs, int nDaughters,
                                            AABBSubdividerConfig const &config)
{
  SplitBVHCostEvaluation result;
  const double minimumDistanceReduction = static_cast<double>(MinimumBVHCostReduction(config));
  if (minimumDistanceReduction < 0.) return result;

  const auto buildConfig = IndexedBVH<BVHBuildReal_t>::VolumeNavigationBuildConfig(nDaughters);
  result.originalBVH =
      std::make_unique<IndexedBVH<BVHBuildReal_t>>(splitBVH.GetRootId(), originalAABBs, nDaughters, buildConfig);
  const double originalCost = EstimatedBVHCost(*result.originalBVH, config);
  const double splitCost    = EstimatedBVHCost(splitBVH, config);
  result.distanceAccepted = originalCost > 0. && (originalCost - splitCost) / originalCost >= minimumDistanceReduction;
  return result;
}

void RejectSplitProduct(BVHBoxProduct &product)
{
  const std::size_t rejectedSplits = product.splitDecisions.size();
  const std::size_t rejectedExtraBoxes =
      product.boxToDaughter.size() > static_cast<std::size_t>(product.nDaughters)
          ? product.boxToDaughter.size() - static_cast<std::size_t>(product.nDaughters)
          : 0;
  product.stats.costRejectedLogicalVolumes += 1;
  product.stats.costRejectedSplits += rejectedSplits;
  product.stats.costRejectedExtraBoxes += rejectedExtraBoxes;
  assert(product.stats.acceptedSplits >= rejectedSplits);
  assert(product.stats.splitBoxes >= rejectedSplits + rejectedExtraBoxes);
  assert(product.stats.extraBoxes >= rejectedExtraBoxes);
  product.stats.acceptedSplits -= rejectedSplits;
  product.stats.splitBoxes -= rejectedSplits + rejectedExtraBoxes;
  product.stats.extraBoxes -= rejectedExtraBoxes;
  product.hasSubdivision = false;
  std::vector<Vector3D<Precision>>{}.swap(product.boxes);
  std::vector<int>{}.swap(product.boxToDaughter);
  std::vector<bvhrecipe::Daughter>{}.swap(product.splitDecisions);
}

void ConstructIndexedBVHRange(std::size_t begin, std::size_t end, void *context)
{
  auto &buildContext = *static_cast<IndexedBVHConstructionContext *>(context);
  for (std::size_t jobIndex = begin; jobIndex < end; ++jobIndex) {
    const std::size_t productIndex = (*buildContext.jobs)[jobIndex];
    auto &boxProduct               = (*buildContext.boxProducts)[productIndex];
    auto &constructionProduct      = (*buildContext.constructionProducts)[productIndex];

    auto start = Clock::now();
    std::unique_ptr<IndexedBVH<BVHBuildReal_t>> bvh;
    if (boxProduct.hasSubdivision) {
      // Build the mapped split hierarchy first: the acceptance model must measure the actual
      // IndexedBVH topology, not infer its cost from the input box count. A Distance rejection
      // restores the original identity hierarchy; an accepted split serves both Distance and Safety.
      const int primitiveCount = static_cast<int>(boxProduct.boxToDaughter.size());
      const auto buildConfig   = IndexedBVH<BVHBuildReal_t>::VolumeNavigationBuildConfig(primitiveCount);
      bvh = std::make_unique<IndexedBVH<BVHBuildReal_t>>(boxProduct.logicalVolume->id(), boxProduct.boxes.data(),
                                                         primitiveCount, boxProduct.boxToDaughter.data(),
                                                         boxProduct.nDaughters, buildConfig);
      if (!boxProduct.fromRecipe && MinimumBVHCostReduction(buildContext.subdivisionConfig) >= 0.) {
        ++boxProduct.stats.costEvaluatedLogicalVolumes;
        auto evaluation =
            EvaluateSplitBVHCost(*bvh, boxProduct.originalAABBs, boxProduct.nDaughters, buildContext.subdivisionConfig);
        if (!evaluation.distanceAccepted) {
          ++boxProduct.stats.distanceCostRejectedLogicalVolumes;
          RejectSplitProduct(boxProduct);
          bvh = std::move(evaluation.originalBVH);
        }
      }
    } else {
      const auto buildConfig = IndexedBVH<BVHBuildReal_t>::VolumeNavigationBuildConfig(boxProduct.nDaughters);
      bvh = std::make_unique<IndexedBVH<BVHBuildReal_t>>(boxProduct.logicalVolume->id(), boxProduct.originalAABBs,
                                                         boxProduct.nDaughters, buildConfig);
    }

    constructionProduct.elapsed        = ElapsedSeconds(start);
    constructionProduct.allocatedBytes = bvh->GetAllocatedSize();
    constructionProduct.bvh            = std::move(bvh);
  }
}

std::vector<IndexedBVHConstructionProduct> ConstructIndexedBVHs(std::vector<BVHBoxProduct> &boxProducts,
                                                                AABBSubdividerConfig const &subdivisionConfig)
{
  std::vector<std::size_t> jobs;
  jobs.reserve(boxProducts.size());
  for (std::size_t index = 0; index < boxProducts.size(); ++index) {
    if (boxProducts[index].nDaughters > 0) jobs.push_back(index);
  }

  std::vector<IndexedBVHConstructionProduct> constructionProducts(boxProducts.size());
  IndexedBVHConstructionContext context{&boxProducts, &jobs, &constructionProducts, subdivisionConfig};
  CurrentBuildExecutor().ParallelFor(0, jobs.size(), 1, &ConstructIndexedBVHRange, &context);
  return constructionProducts;
}

void CollectIndexedBVHProduct(BVHBoxProduct const &boxProduct, IndexedBVHConstructionProduct &constructionProduct,
                              IndexedBVHOwnerTable &nextBVHs, BVHBuildStatistics &stats)
{
  AccumulateProductStats(stats, boxProduct.stats);

  if (boxProduct.nDaughters == 0) return;

  stats.finalBoxes += boxProduct.hasSubdivision ? boxProduct.boxToDaughter.size() : boxProduct.nDaughters;
  if (boxProduct.hasSubdivision) ++stats.splitLogicalVolumes;
  stats.indexedBVHConstructionWorkerElapsed += constructionProduct.elapsed;
  stats.indexedHostAllocatedBytes += constructionProduct.allocatedBytes;
  stats.bvhConstructionWorkerElapsed += constructionProduct.elapsed;
  stats.hostAllocatedBytes += constructionProduct.allocatedBytes;
  assert(constructionProduct.bvh);
  nextBVHs[boxProduct.logicalVolume->id()] = std::move(constructionProduct.bvh);
}

void ReplacePublishedIndexedBVHs(IndexedBVHOwnerTable nextBVHs)
{
  std::vector<IndexedBVH<BVHBuildReal_t> *> nextLookup(nextBVHs.size(), nullptr);
  for (std::size_t index = 0; index < nextBVHs.size(); ++index)
    nextLookup[index] = nextBVHs[index].get();

  gOwnedIndexedBVHs.swap(nextBVHs);
  hBVH<BVHBuildReal_t>.swap(nextLookup);
}

bool CreateCurrentRecipe(bvhrecipe::Recipe &recipe, std::string &error)
{
  auto const *world = GeoManager::Instance().GetWorld();
  if (!gHasSuccessfulBuild || !world || world != gLastBuildWorld) {
    error = "the current geometry has not completed BVHManager::Init()";
    return false;
  }

  BVHBuildStatistics ignoredStats;
  auto const inputs     = SnapshotBuildInputs(gLastBuildLogicalVolumes, ignoredStats);
  auto const identities = ComputeRecipeIdentities(inputs);

  bvhrecipe::Recipe built;
  built.vecgeomVersion      = vecgeom_version;
  built.precisionBytes      = sizeof(Precision);
  built.bvhRealBytes        = sizeof(BVHBuildReal_t);
  built.worldLabel          = world->GetLabel();
  built.geometryFingerprint = identities.geometryFingerprint;

  for (std::size_t inputIndex = 0; inputIndex < inputs.size(); ++inputIndex) {
    auto const &input = inputs[inputIndex];
    if (input.nDaughters == 0) continue;
    if (input.logicalVolume->id() >= hBVH<BVHBuildReal_t>.size() || !hBVH<BVHBuildReal_t>[input.logicalVolume->id()]) {
      error = "missing BVH for logical volume '" + input.logicalVolume->GetLabel() + "'";
      return false;
    }
    if (input.logicalVolume->id() > static_cast<unsigned>(std::numeric_limits<int>::max())) {
      error = "logical volume id exceeds the recipe format limit";
      return false;
    }

    bvhrecipe::Volume volume;
    volume.logicalId     = static_cast<int>(input.logicalVolume->id());
    volume.label         = input.logicalVolume->GetLabel();
    volume.fingerprint   = identities.volumes[inputIndex].fingerprint;
    volume.daughterCount = input.nDaughters;
    auto const decisions = gLastBuildSplitDecisions.find(input.logicalVolume);
    if (decisions != gLastBuildSplitDecisions.end()) volume.splitDaughters = decisions->second;
    built.volumes.push_back(std::move(volume));
  }

  recipe = std::move(built);
  error.clear();
  return true;
}

void PublishMissingEnvironmentRecipe()
{
  if (gSubdivisionRecipeCreateFile.empty()) return;

  // Clear the cached environment path so a later Init() reads the file created
  // here, or retries creation after a publication failure.
  std::string const fileName = std::move(gSubdivisionRecipeCreateFile);
  gSubdivisionRecipeCreateFile.clear();
  gSubdivisionRecipeEnvironmentFile.clear();

  auto const start = Clock::now();
  bvhrecipe::Recipe recipe;
  std::string error;
  bvhrecipe::PublishResult result = bvhrecipe::PublishResult::kError;
  if (CreateCurrentRecipe(recipe, error)) result = bvhrecipe::WriteIfAbsent(fileName, recipe, error);
  gBuildStats.recipeAutoWriteElapsed = ElapsedSeconds(start);

  switch (result) {
  case bvhrecipe::PublishResult::kCreated:
    ++gBuildStats.recipeAutoWriteCreated;
    VECGEOM_LOG(info) << "Created BVH AABB subdivision recipe '" << fileName << "'";
    break;
  case bvhrecipe::PublishResult::kExists:
    ++gBuildStats.recipeAutoWriteSkipped;
    VECGEOM_LOG(info) << "Did not replace BVH AABB subdivision recipe '" << fileName
                      << "' created by another initializer";
    break;
  case bvhrecipe::PublishResult::kError:
    ++gBuildStats.recipeAutoWriteFailed;
    VECGEOM_LOG(warning) << "Cannot create BVH AABB subdivision recipe '" << fileName << "': " << error;
    break;
  }
}

} // namespace
#endif

void BVHManager::SetAABBSubdivisionMode(AABBSubdivisionMode mode)
{
#ifndef VECCORE_CUDA
  gSubdivisionConfig.mode  = mode;
  gSubdivisionModeExplicit = true;
#else
  (void)mode;
#endif
}

AABBSubdivisionMode BVHManager::GetAABBSubdivisionMode()
{
#ifndef VECCORE_CUDA
  LoadSubdivisionEnvironment();
  return gSubdivisionConfig.mode;
#else
  return AABBSubdivisionMode::kOff;
#endif
}

void BVHManager::SetAABBSubdivisionMaxM(int maxM)
{
#ifndef VECCORE_CUDA
  gSubdivisionMaxM         = maxM > 0 ? std::min(maxM, kMaximumAABBSubdivisionM) : 0;
  gSubdivisionMaxMExplicit = true;
#else
  (void)maxM;
#endif
}

int BVHManager::GetAABBSubdivisionMaxM()
{
#ifndef VECCORE_CUDA
  LoadSubdivisionEnvironment();
  return gSubdivisionMaxM > 0 ? std::min(gSubdivisionMaxM, kMaximumAABBSubdivisionM) : 0;
#else
  return 0;
#endif
}

void BVHManager::SetAABBSubdivisionConfig(AABBSubdividerConfig const &config)
{
#ifndef VECCORE_CUDA
  gSubdivisionConfig       = config;
  gSubdivisionModeExplicit = true;
  gSubdivisionMaxM         = 0;
  gSubdivisionMaxMExplicit = true;
#else
  (void)config;
#endif
}

AABBSubdividerConfig BVHManager::GetAABBSubdivisionConfig()
{
#ifndef VECCORE_CUDA
  return CurrentSubdivisionConfig();
#else
  AABBSubdividerConfig config;
  config.mode = AABBSubdivisionMode::kOff;
  return config;
#endif
}

BVHBuildStatistics BVHManager::GetBuildStatistics()
{
#ifndef VECCORE_CUDA
  return gBuildStats;
#else
  return {};
#endif
}

bool BVHManager::ReadAABBSubdivisionRecipe(std::string const &fileName, std::string *errorMessage)
{
#ifndef VECCORE_CUDA
  gSubdivisionRecipeExplicit = true;
  gSubdivisionRecipeEnvironmentFile.clear();
  gSubdivisionRecipeCreateFile.clear();
  gLoadedRecipe.reset();

  auto recipe = std::make_unique<bvhrecipe::Recipe>();
  std::string error;
  bool const success = bvhrecipe::Read(fileName, *recipe, error);
  if (success) gLoadedRecipe = std::move(recipe);
  if (errorMessage) *errorMessage = error;
  return success;
#else
  (void)fileName;
  if (errorMessage) *errorMessage = "AABB subdivision recipes are a host-only build feature";
  return false;
#endif
}

bool BVHManager::WriteAABBSubdivisionRecipe(std::string const &fileName, std::string *errorMessage)
{
#ifndef VECCORE_CUDA
  bvhrecipe::Recipe recipe;
  std::string error;
  bool success = CreateCurrentRecipe(recipe, error);
  if (success) success = bvhrecipe::Write(fileName, recipe, error);
  if (errorMessage) *errorMessage = error;
  return success;
#else
  (void)fileName;
  if (errorMessage) *errorMessage = "AABB subdivision recipes are a host-only build feature";
  return false;
#endif
}

void BVHManager::ClearAABBSubdivisionRecipe()
{
#ifndef VECCORE_CUDA
  gSubdivisionRecipeExplicit = true;
  gSubdivisionRecipeEnvironmentFile.clear();
  gSubdivisionRecipeCreateFile.clear();
  gLoadedRecipe.reset();
#endif
}

void BVHManager::SetBuildExecutor(BVHBuildExecutor *executor)
{
#ifndef VECCORE_CUDA
  gBuildExecutor = executor;
#else
  (void)executor;
#endif
}

void BVHManager::SetBuildThreadCount(unsigned nThreads)
{
#ifndef VECCORE_CUDA
  gBuildThreadCount = nThreads;
#else
  (void)nThreads;
#endif
}

unsigned BVHManager::GetBuildThreadCount()
{
#ifndef VECCORE_CUDA
  return gBuildThreadCount;
#else
  return 1;
#endif
}

void BVHManager::Init()
{
#ifndef VECCORE_CUDA
  auto initStart = Clock::now();
  BVHBuildStatistics nextBuildStats;
  LoadSubdivisionRecipeEnvironment();
#endif

  std::vector<LogicalVolume const *> lvols;
  GeoManager::Instance().GetAllLogicalVolumes(lvols);
#ifndef VECCORE_CUDA
  auto subdivisionConfig               = CurrentSubdivisionConfig();
  auto inputs                          = SnapshotBuildInputs(lvols, nextBuildStats);
  auto boxPreparationStart             = Clock::now();
  auto products                        = BuildBoxProducts(inputs, subdivisionConfig, nextBuildStats);
  nextBuildStats.boxPreparationElapsed = ElapsedSeconds(boxPreparationStart);

  auto indexedConstructionStart                = Clock::now();
  auto indexedConstructionProducts             = ConstructIndexedBVHs(products, subdivisionConfig);
  nextBuildStats.indexedBVHConstructionElapsed = ElapsedSeconds(indexedConstructionStart);
  nextBuildStats.bvhConstructionElapsed        = nextBuildStats.indexedBVHConstructionElapsed;
  // There may be volumes not used in the hierarchy, so the maximum index may be larger.
  IndexedBVHOwnerTable nextIndexedBVHs(GeoManager::Instance().GetLogicalVolumesMap().size());
  for (std::size_t index = 0; index < products.size(); ++index) {
    CollectIndexedBVHProduct(products[index], indexedConstructionProducts[index], nextIndexedBVHs, nextBuildStats);
  }

  decltype(gLastBuildSplitDecisions) nextSplitDecisions;
  nextSplitDecisions.reserve(std::count_if(
      products.begin(), products.end(), [](BVHBoxProduct const &product) { return !product.splitDecisions.empty(); }));
  for (auto &product : products) {
    if (!product.splitDecisions.empty())
      nextSplitDecisions.emplace(product.logicalVolume, std::move(product.splitDecisions));
  }
  nextBuildStats.recipeDecisionBytesRetained =
      nextSplitDecisions.bucket_count() * sizeof(void *) +
      nextSplitDecisions.size() * sizeof(decltype(nextSplitDecisions)::value_type);
  for (auto const &entry : nextSplitDecisions)
    nextBuildStats.recipeDecisionBytesRetained += entry.second.capacity() * sizeof(bvhrecipe::Daughter);

  auto nextLogicalVolumes = lvols;
  auto *nextWorld         = GeoManager::Instance().GetWorld();
  ReplacePublishedIndexedBVHs(std::move(nextIndexedBVHs));
  gLastBuildSplitDecisions.swap(nextSplitDecisions);
  gLastBuildLogicalVolumes.swap(nextLogicalVolumes);
  gLastBuildWorld     = nextWorld;
  gHasSuccessfulBuild = true;
  gBuildStats         = nextBuildStats;
  PublishMissingEnvironmentRecipe();
#else
  // There may be volumes not used in the hierarchy, so the maximum index may be larger.
  hBVH<Real_t>.resize(GeoManager::Instance().GetLogicalVolumesMap().size());
  for (auto logical_volume : lvols) {
    int n{0};
    auto ptrAABB = ABBoxManager<Precision>::Instance().GetABBoxes(logical_volume, n);
    if (n == 0) {
      hBVH<Real_t>[logical_volume->id()] = nullptr;
      continue;
    }
    hBVH<Real_t>[logical_volume->id()] =
        new IndexedBVH<Real_t>(logical_volume->id(), ptrAABB, n, IndexedBVH<Real_t>::VolumeNavigationBuildConfig(n));
  }
#endif

#ifndef VECCORE_CUDA
  gBuildStats.initElapsed = ElapsedSeconds(initStart);
#endif
}

cuda::IndexedBVH<BVHManager::Real_t> const *BVHManager::DeviceInit()
{
#ifdef VECGEOM_CUDA_INTERFACE
#ifndef VECCORE_CUDA
  auto indexedDeviceStart                 = Clock::now();
  gBuildStats.deviceAllocatedBytes        = 0;
  gBuildStats.deviceAllocationCount       = 0;
  gBuildStats.indexedDeviceAllocatedBytes = 0;
#endif

  using HostBVH = IndexedBVH<Real_t>;
  using Layout  = typename HostBVH::DeviceImageLayout;
  static_assert(sizeof(Layout) <= sizeof(HostBVH), "IndexedBVH device layout must fit in one object slot");

  const std::size_t n = hBVH<Real_t>.size();
  if (n == 0) return nullptr;

  // Reserve every array in one aligned arena. Object slots initially hold descriptors; the
  // initialization kernel consumes and replaces them with non-owning IndexedBVH views.
  std::vector<Layout> layouts(n);
  std::size_t arenaSize = n * sizeof(HostBVH);
  for (std::size_t id = 0; id < n; ++id) {
    auto const *hostBVH = hBVH<Real_t>[id];
    if (!hostBVH) continue;
    layouts[id] = hostBVH->MakeDeviceImageLayout(arenaSize);
  }

  std::vector<char> hostArena(arenaSize, 0);
  for (std::size_t id = 0; id < n; ++id) {
    auto const *hostBVH = hBVH<Real_t>[id];
    if (!hostBVH) continue;
    hostBVH->PackDeviceImage(hostArena.data(), hostArena.size(), id * sizeof(HostBVH), layouts[id]);
  }

  auto deviceArena = vecgeom::cxx::CudaManager::Instance().AllocateDeviceBuffer(arenaSize);
  deviceArena.ToDevice(hostArena.data(), arenaSize);
  auto *ptr =
      vecgeom::cuda::InitializeDeviceIndexedBVHBuffer<Real_t>(deviceArena.GetPtr(), n, sizeof(HostBVH), sizeof(Layout));

#ifndef VECCORE_CUDA
  gBuildStats.deviceAllocatedBytes        = arenaSize;
  gBuildStats.deviceAllocationCount       = 1;
  gBuildStats.indexedDeviceAllocatedBytes = arenaSize;
  gBuildStats.indexedDeviceInitElapsed    = ElapsedSeconds(indexedDeviceStart);
  gBuildStats.deviceInitElapsed           = gBuildStats.indexedDeviceInitElapsed;
#endif

  return ptr;
#else
  VECGEOM_LOG(error) << "Cannot initialize BVH device: CUDA is not configured";
  return nullptr;
#endif
}

cuda::IndexedBVH<BVHManager::Real_t> const *BVHManager::GetDeviceBVH()
{
#ifdef VECGEOM_CUDA_INTERFACE
  return cuda::GetDeviceIndexedBVHBuffer<Real_t>();
#else
  return nullptr;
#endif
}

} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom
