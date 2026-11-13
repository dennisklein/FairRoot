/// \file BVHNavigatorV.h
/// \brief IndexedBVH-accelerated navigator behind the virtual VNavigator interface.

#ifndef VECGEOM_NAVIGATION_BVHNAVIGATORV_H_
#define VECGEOM_NAVIGATION_BVHNAVIGATORV_H_

#include "VecGeom/base/IndexedBVH.h"
#include "VecGeom/management/BVHManager.h"
#include "VecGeom/navigation/BVHNavigator.h"
#include "VecGeom/navigation/BVHSafetyEstimator.h"
#include "VecGeom/navigation/VNavigator.h"
#include "VecGeom/volumes/LogicalVolume.h"
#include "VecGeom/volumes/PlacedVolume.h"

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

/**
 * @brief IndexedBVH navigator behind the virtual VNavigator interface.
 *
 * The VNavigator counterpart of the standalone BVHNavigator, for consumers that assign a
 * navigator per logical volume. Pair it with BVHLevelLocator and BVHSafetyEstimator.
 * Only the daughter query is implemented here; VNavigatorHelper provides the rest.
 * BVHManager::Init() must have run after GeoManager::CloseGeometry(), or the hierarchies do not exist.
 */
template <bool MotherIsConvex = false>
class BVHNavigatorV : public VNavigatorHelper<BVHNavigatorV<MotherIsConvex>, MotherIsConvex> {

private:
  /** Constructor. Private since this is a singleton class accessed only via Instance(). */
  BVHNavigatorV() : VNavigatorHelper<BVHNavigatorV<MotherIsConvex>, MotherIsConvex>(BVHSafetyEstimator::Instance()) {}

public:
  using Base = VNavigatorHelper<BVHNavigatorV<MotherIsConvex>, MotherIsConvex>;
  using Base::PrepareOutState;

  /// Safety from the same hierarchy.
  typedef BVHSafetyEstimator SafetyEstimator_t;

  static constexpr const char *gClassNameString = "BVHNavigatorV";

  /** Returns the instance of this singleton class. */
  static VNavigator *Instance()
  {
    static BVHNavigatorV instance;
    return &instance;
  }

  /**
   * Nearest daughter of @p lvol hit within @p step. The volume last exited in @p in_state is refused only at a
   * non-positive distance, so a non-convex daughter the ray leaves can be re-entered further on.
   * @returns Always false: the outgoing state is not modified.
   */
  VECCORE_ATT_HOST_DEVICE
  bool CheckDaughterIntersections(LogicalVolume const *lvol, Vector3D<Precision> const &localpoint,
                                  Vector3D<Precision> const &localdir, NavigationState const *in_state,
                                  NavigationState * /* out_state */, Precision &step,
                                  VPlacedVolume const *&hitcandidate) const final
  {
    auto const *bvh = BVHManager::GetBVH(lvol->id());
    if (bvh == nullptr) return false;

    // A placed-volume id is unsigned; -1 as a long is the "no exclusion" marker.
    long last_exited_id          = -1;
    VPlacedVolume const *lastvol = in_state ? in_state->GetLastExited() : nullptr;
    if (lastvol != nullptr) last_exited_id = lastvol->id();

    long hitcandidate_index = -1;
    BVHNavigator::CheckDaughterIntersectionsBlocked<BVHNavigator>(*bvh, localpoint, localdir, step, last_exited_id,
                                                                  hitcandidate_index);
    if (hitcandidate_index >= 0) hitcandidate = (*lvol->GetDaughtersp())[hitcandidate_index];

    return false; /* return value indicates whether out_state has been modified */
  }

  /**
   * Nearest daughter of @p lvol hit within @p step. At a non-positive distance it refuses @p blocked and any
   * daughter the ray is leaving, as SimpleABBoxNavigator does. This is the query of Geant4 glue.
   * @returns Always false: no outgoing state is modified.
   */
  VECCORE_ATT_HOST_DEVICE
  bool CheckDaughterIntersections(LogicalVolume const *lvol, Vector3D<Precision> const &localpoint,
                                  Vector3D<Precision> const &localdir, VPlacedVolume const *blocked, Precision &step,
                                  VPlacedVolume const *&hitcandidate) const override
  {
    auto const *bvh = BVHManager::GetBVH(lvol->id());
    if (bvh == nullptr) return false;

    // A placed-volume id is unsigned; -1 as a long is the "no exclusion" marker.
    long last_exited_id = (blocked != nullptr) ? static_cast<long>(blocked->id()) : -1;

    long hitcandidate_index = -1;
    BVHNavigator::CheckDaughterIntersectionsBlocked<BVHNavigator>(*bvh, localpoint, localdir, step, last_exited_id,
                                                                  hitcandidate_index);
    if (hitcandidate_index >= 0) hitcandidate = (*lvol->GetDaughtersp())[hitcandidate_index];

    return false; /* return value indicates whether out_state has been modified */
  }
}; // end class declaration

} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif // VECGEOM_NAVIGATION_BVHNAVIGATORV_H_
