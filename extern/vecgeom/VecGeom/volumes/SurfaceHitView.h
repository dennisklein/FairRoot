/// @file SurfaceHitView.h
/// @brief Optional exit-normal information returned by distance queries.

#ifndef VECGEOM_VOLUMES_SURFACEHITVIEW_H_
#define VECGEOM_VOLUMES_SURFACEHITVIEW_H_

#include "VecGeom/base/Global.h"
#include "VecGeom/base/Vector3D.h"

namespace vecgeom {

inline namespace VECGEOM_IMPL_NAMESPACE {

/// @brief Optional exit information associated with a scalar `DistanceToOut` result.
/// @details Normal storage is supplied by the caller. When `fLeftBehind` is
/// true, the complete solid lies behind or on the plane through the hit point
/// with the returned outward normal. The normal is usable only when
/// `fLeftBehind` is true. This flag is not the validity result returned by
/// VecGeom's standalone `Normal()` query.
template <typename Real_t = Precision>
struct SurfaceHitView {
  Vector3D<Real_t> *fNormal = nullptr; ///< Optional destination for a usable unit outward normal.
  bool fLeftBehind          = false;   ///< Whether the hit plane leaves the complete solid behind or on it.

  /// @brief Check whether normal output was requested.
  /// @return True when `fNormal` points to writable storage.
  VECCORE_ATT_HOST_DEVICE
  bool WantsNormal() const { return fNormal != nullptr; }

  /// @brief Reset validity without touching caller-owned normal storage.
  VECCORE_ATT_HOST_DEVICE
  void Clear() { fLeftBehind = false; }

  /// @brief Store whether the returned normal defines a supporting exit plane.
  /// @param leftBehind True when the complete solid lies behind or on the hit plane.
  VECCORE_ATT_HOST_DEVICE
  void SetLeftBehind(bool leftBehind) { fLeftBehind = leftBehind; }

  /// @brief Store a unit outward normal if requested.
  /// @param normal Unit outward normal in the query-local coordinate frame.
  VECCORE_ATT_HOST_DEVICE
  void SetNormal(Vector3D<Real_t> const &normal) const
  {
    if (fNormal) *fNormal = normal;
  }
};

} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif // VECGEOM_VOLUMES_SURFACEHITVIEW_H_
