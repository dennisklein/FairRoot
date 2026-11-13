#ifndef VECGEOM_VOLUMES_BOOLEANSTRUCT_H_
#define VECGEOM_VOLUMES_BOOLEANSTRUCT_H_
#include "VecGeom/base/Global.h"
#include <cstdint>

namespace vecgeom {

// Declare types shared by cxx and cuda.
enum BooleanOperation { kUnion, kIntersection, kSubtraction };

inline namespace VECGEOM_IMPL_NAMESPACE {

/**
 * A class representing a simple UNPLACED boolean volume A-B
 * It takes two template arguments:
 * 1.: the mother (or left) volume A in unplaced form
 * 2.: the (or right) volume B in placed form, acting on A with a boolean operation;
 * the placement is with respect to the left volume
 */
struct BooleanStruct {
  VPlacedVolume const *fLeftVolume;
  VPlacedVolume const *fRightVolume;
  BooleanOperation const fOp;
  std::uint16_t fLeftScore       = 1; ///< Cached left Inside cost in approximate box units, saturated at 65535.
  std::uint16_t fRightScore      = 1; ///< Cached right Inside cost in approximate box units, saturated at 65535.
  std::uint8_t fIdentityMask     = 0; ///< Exact identity placements: bit 0 left, bit 1 right; zero disables shortcuts.
  mutable Precision fCapacity    = -1;
  mutable Precision fSurfaceArea = -1;

  VECCORE_ATT_HOST_DEVICE
  BooleanStruct(BooleanOperation op, VPlacedVolume const *left, VPlacedVolume const *right)
      : fLeftVolume(left), fRightVolume(right), fOp(op)
  {
#ifndef VECCORE_CUDA
    UpdateScores();
#endif
  }

#ifndef VECCORE_CUDA
  /// @brief Refresh operand costs and placement flags after construction or a subtree rewrite.
  /// @details Children must already be updated. Adds cached child costs without traversing the Boolean tree.
  /// Leaves with existing bounds checks retain unit cost; polygon scans and polycone lookup get coarse weights.
  void UpdateScores();
#endif

  /// @brief Refresh exact-identity flags after operand replacement; placements must be finalized.
  VECCORE_ATT_HOST_DEVICE
  void UpdatePlacementFlags();
}; // End struct

} // namespace VECGEOM_IMPL_NAMESPACE

} // namespace vecgeom

#endif /* VECGEOM_VOLUMES_BOOLEANSTRUCT_H_ */
