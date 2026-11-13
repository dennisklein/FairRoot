//===-- source/BVHSubdivisionRecipe.h -----------------------------*- C++ -*-===//
//
/// @file BVHSubdivisionRecipe.h
/// @brief Private representation for persisted BVH slicing decisions.
//
//===----------------------------------------------------------------------===//

#ifndef VECGEOM_SOURCE_BVHSUBDIVISIONRECIPE_H_
#define VECGEOM_SOURCE_BVHSUBDIVISIONRECIPE_H_

#include "VecGeom/base/Global.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {
namespace bvhrecipe {

/// Increment when the fixed little-endian serialized grammar changes incompatibly.
constexpr int kFormatVersion = 1;

/// Fixed uniform split decision for one daughter.
struct Daughter {
  int daughterId{-1};
  int axis{-1}; ///< Mother-frame split axis: 0=X, 1=Y, 2=Z.
  int M{1};     ///< Number of uniform slabs along `axis`.
};

/// Complete set of split decisions for one logical volume.
///
/// Daughters absent from `splitDaughters` intentionally retain their original
/// AABB. Therefore a matched volume bypasses live heuristic evaluation even
/// when this list is empty.
struct Volume {
  int logicalId{-1};
  std::string label;
  std::uint64_t fingerprint{0};
  std::size_t daughterCount{0};
  std::vector<Daughter> splitDaughters;
};

/// Versioned, geometry-bound set of fixed mother BVH split decisions.
struct Recipe {
  std::string vecgeomVersion;           ///< Diagnostic producer version; format version gates grammar.
  int precisionBytes{0};                ///< Host geometry precision width.
  int bvhRealBytes{0};                  ///< BVH coordinate precision width.
  std::string worldLabel;               ///< Top-level geometry identity.
  std::uint64_t geometryFingerprint{0}; ///< Complete reachable-geometry fingerprint.
  std::vector<Volume> volumes;          ///< One complete decision record per mother.
};

/// Result of publishing a recipe without replacing an existing recipe file.
enum class PublishResult {
  kCreated, ///< This call created the requested recipe file.
  kExists,  ///< Another producer or an earlier call already created it.
  kError    ///< Publication failed; `error` describes the failure.
};

/// Parse and structurally validate one complete recipe.
///
/// @param fileName Recipe file path to read.
/// @param[out] recipe Parsed product, replaced only after complete success.
/// @param[out] error Parse or I/O diagnostic, cleared on success.
/// @return True when the complete supported recipe file was parsed.
bool Read(std::string const &fileName, Recipe &recipe, std::string &error);

/// Write one complete compact split-decision recipe.
///
/// Integer fields are emitted byte-by-byte in little-endian order; no C++
/// object representation or host byte order is serialized.
///
/// @param fileName Recipe file path to create or replace.
/// @param recipe Complete product to serialize.
/// @param[out] error I/O diagnostic, cleared on success.
/// @return True when all data was flushed and the file closed successfully.
bool Write(std::string const &fileName, Recipe const &recipe, std::string &error);

/// Publish one complete recipe only when the destination is absent.
///
/// The payload is first written beside the destination, then published without
/// replacing an existing file. This supports environment-controlled cache
/// creation by concurrent initializers without silently clobbering a winner.
/// Filesystems without atomic hard-link publication return `kError`.
///
/// @param fileName Recipe file path to create if absent.
/// @param recipe Complete product to serialize.
/// @param[out] error Publication diagnostic, cleared for `kCreated` and `kExists`.
/// @return Whether this call created the file, observed an existing file, or failed.
PublishResult WriteIfAbsent(std::string const &fileName, Recipe const &recipe, std::string &error);

} // namespace bvhrecipe
} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif // VECGEOM_SOURCE_BVHSUBDIVISIONRECIPE_H_
