// This file is part of VecGeom and is distributed under the
// conditions in the file LICENSE.txt in the top directory.

/// @file BooleanFactory.h
/// @brief Build equivalent Boolean solids with accelerated MultiUnion components.

#ifndef VECGEOM_MANAGEMENT_BOOLEANFACTORY_H_
#define VECGEOM_MANAGEMENT_BOOLEANFACTORY_H_

#include "VecGeom/base/Global.h"
#include <unordered_map>

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

class LogicalVolume;
class Transformation3D;
class VPlacedVolume;
class VUnplacedVolume;

/**
 * @brief Convert sufficiently large union groups and subtraction cutters to MultiUnions.
 * @details Input solids are unchanged. Repeated inputs share their replacement within
 * this factory. Created shapes, logical volumes, placements and transforms follow the
 * usual geometry lifetime; destroying the factory only releases its lookup tables.
 * Results can reuse input geometry; keep that geometry alive and unchanged
 * while using the factory or its results.
 * External adapters should call Create on the complete solid before publishing
 * their final logical volume or placement. This factory does not update existing
 * wrappers or make geometry compaction safe for externally retained placements.
 */
class BooleanFactory {
  size_t fMinComponents;
  std::unordered_map<VUnplacedVolume const *, VUnplacedVolume const *> fResults;
  std::unordered_map<VUnplacedVolume const *, LogicalVolume *> fLogicalVolumes;

  VPlacedVolume const *Place(VUnplacedVolume const *solid, Transformation3D const &transform);
  VPlacedVolume const *ConvertPlaced(VPlacedVolume const *placed);

public:
  /**
   * @brief Configure the minimum number of resulting MultiUnion components.
   * @param minComponents Minimum component count; zero disables conversion.
   * @details A MultiUnion always needs at least two components. This threshold
   * counts components, independently of the operand Inside-cost scores.
   */
  explicit BooleanFactory(size_t minComponents) : fMinComponents(minComponents) {}

  /**
   * @brief Build an equivalent solid, reusing previous results for shared inputs.
   * @param solid Input solid, or null.
   * @return Replacement solid, the original pointer if unchanged, or null for null input.
   * @details Flatten union groups and left subtraction spines only when the
   * resulting union has enough finite components. Right subtraction and
   * intersection grouping are retained; their children can be converted.
   * Existing MultiUnions and scaled shapes retain their outer operation.
   */
  VUnplacedVolume const *Create(VUnplacedVolume const *solid);
};

} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif
