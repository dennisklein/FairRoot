/// \file BVHLevelLocator.h
/// \brief Level locator accelerated by the IndexedBVH of each logical volume.

#ifndef VECGEOM_NAVIGATION_BVHLEVELLOCATOR_H_
#define VECGEOM_NAVIGATION_BVHLEVELLOCATOR_H_

#include "VecGeom/base/IndexedBVH.h"
#include "VecGeom/management/BVHManager.h"
#include "VecGeom/navigation/SimpleLevelLocator.h"
#include "VecGeom/navigation/VLevelLocator.h"
#include "VecGeom/volumes/LogicalVolume.h"
#include "VecGeom/volumes/PlacedVolume.h"

#include <string>

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

/**
 * @brief Level locator enumerating daughter candidates through the IndexedBVH of the mother.
 *
 * Candidates are accepted by the CheckCandidateVol() helpers of the other level locators. The
 * first accepted candidate wins, so on a face shared by two daughters the answer depends on the
 * traversal order and may differ from TSimpleLevelLocator.
 */
template <bool IsAssemblyAware = false>
class TBVHLevelLocator : public VLevelLocator {

private:
  /** Constructor. Private since this is a singleton class accessed only via GetInstance(). */
  TBVHLevelLocator() = default;

  /**
   * First daughter of @p lvol whose box contains @p localpoint and which the candidate helpers accept.
   * A daughter split over several boxes may be tested more than once, with the same answer.
   */
  template <bool ExclV, bool ModifyState, bool WithDirection, typename Real_b>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE bool LevelLocateKernel(
      IndexedBVH<Real_b> const &bvh, LogicalVolume const *lvol, VPlacedVolume const *exclvol,
      Vector3D<Precision> const &localpoint, Vector3D<Precision> const *localdir, NavigationState *state,
      VPlacedVolume const *&pvol, Vector3D<Precision> &daughterlocalpoint) const
  {
    const Vector3D<Real_b> point(static_cast<Real_b>(localpoint[0]), static_cast<Real_b>(localpoint[1]),
                                 static_cast<Real_b>(localpoint[2]));
    if (!bvh.NodeContains(0, point)) return false;

    auto const *nodes     = bvh.GetNodes();
    auto const *primId    = bvh.GetPrimId();
    auto const &daughters = *lvol->GetDaughtersp();
    const int nDaughters  = static_cast<int>(daughters.size());

    unsigned int stack[IndexedBVH<Real_b>::kTraversalStackSize];
    int stackSize    = 0;
    unsigned int top = 0;
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
        const int daughterId = bvh.GetDaughterId(boxId);
        if (daughterId < 0 || daughterId >= nDaughters) continue;
        VPlacedVolume const *daughter = daughters[daughterId];
        if (ExclV && daughter == exclvol) continue;
        if constexpr (WithDirection) {
          if (CheckCandidateVolWithDirection<IsAssemblyAware, ModifyState>(daughter, localpoint, *localdir, state, pvol,
                                                                           daughterlocalpoint))
            return true;
        } else {
          if (CheckCandidateVol<IsAssemblyAware, ModifyState>(daughter, localpoint, state, pvol, daughterlocalpoint))
            return true;
        }
      }

    pop:
      if (stackSize == 0) return false;
      top = stack[--stackSize];
    }
  }

  /// Fetches the hierarchy of @p lvol and runs the kernel; volumes without one locate nothing.
  template <bool ExclV, bool ModifyState, bool WithDirection>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE bool Dispatch(LogicalVolume const *lvol, VPlacedVolume const *exclvol,
                                                             Vector3D<Precision> const &localpoint,
                                                             Vector3D<Precision> const *localdir,
                                                             NavigationState *state, VPlacedVolume const *&pvol,
                                                             Vector3D<Precision> &daughterlocalpoint) const
  {
    auto const *bvh = BVHManager::GetBVH(lvol->id());
    if (bvh == nullptr) return false;
    return LevelLocateKernel<ExclV, ModifyState, WithDirection>(*bvh, lvol, exclvol, localpoint, localdir, state, pvol,
                                                                daughterlocalpoint);
  }

public:
  /** Returns the instance of this singleton class. */
  static VLevelLocator const *GetInstance()
  {
    static TBVHLevelLocator instance;
    return &instance;
  }

  std::string GetName() const final { return IsAssemblyAware ? "AssemblyAwareBVHLevelLocator" : "BVHLevelLocator"; }

  VECCORE_ATT_HOST_DEVICE
  bool LevelLocate(LogicalVolume const *lvol, Vector3D<Precision> const &localpoint, VPlacedVolume const *&pvol,
                   Vector3D<Precision> &daughterlocalpoint) const final
  {
    return Dispatch<false, false, false>(lvol, nullptr, localpoint, nullptr, nullptr, pvol, daughterlocalpoint);
  }

  VECCORE_ATT_HOST_DEVICE
  bool LevelLocate(LogicalVolume const *lvol, Vector3D<Precision> const &localpoint, NavigationState &state,
                   Vector3D<Precision> &daughterlocalpoint) const final
  {
    VPlacedVolume const *pvol = nullptr;
    return Dispatch<false, true, false>(lvol, nullptr, localpoint, nullptr, &state, pvol, daughterlocalpoint);
  }

  VECCORE_ATT_HOST_DEVICE
  bool LevelLocateExclVol(LogicalVolume const *lvol, VPlacedVolume const *exclvol,
                          Vector3D<Precision> const &localpoint, VPlacedVolume const *&pvol,
                          Vector3D<Precision> &daughterlocalpoint) const final
  {
    return Dispatch<true, false, false>(lvol, exclvol, localpoint, nullptr, nullptr, pvol, daughterlocalpoint);
  }

  VECCORE_ATT_HOST_DEVICE
  bool LevelLocateExclVol(LogicalVolume const *lvol, VPlacedVolume const *exclvol,
                          Vector3D<Precision> const &localpoint, Vector3D<Precision> const &localdirection,
                          VPlacedVolume const *&pvol, Vector3D<Precision> &daughterlocalpoint) const final
  {
    return Dispatch<true, false, true>(lvol, exclvol, localpoint, &localdirection, nullptr, pvol, daughterlocalpoint);
  }
}; // end class declaration

using BVHLevelLocator              = TBVHLevelLocator<>;
using BVHAssemblyAwareLevelLocator = TBVHLevelLocator<true>;

} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif // VECGEOM_NAVIGATION_BVHLEVELLOCATOR_H_
