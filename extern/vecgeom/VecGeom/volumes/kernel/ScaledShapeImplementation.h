/// @file ScaledShapeImplementation.h
/// @brief Navigation kernels for a diagonally scaled constituent solid.
/// @author Mihaela Gheata (mihaela.gheata@cern.ch)

#ifndef VECGEOM_VOLUMES_KERNEL_SCALEDSHAPEIMPLEMENTATION_H_
#define VECGEOM_VOLUMES_KERNEL_SCALEDSHAPEIMPLEMENTATION_H_

#include "VecGeom/volumes/kernel/GenericKernels.h"
#include "VecGeom/base/Vector3D.h"
#include "VecGeom/volumes/ScaledShapeStruct.h"
#include "VecGeom/volumes/SurfaceHitView.h"

#include <cstdio>

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE(struct ScaledShapeImplementation;);
VECGEOM_DEVICE_DECLARE_CONV(struct, ScaledShapeImplementation);

inline namespace VECGEOM_IMPL_NAMESPACE {

class PlacedScaledShape;
class UnplacedScaledShape;

// template <typename T>
// struct ScaledShapeStruct;

/// @brief Delegate navigation to a constituent in its unscaled frame.
struct ScaledShapeImplementation {

  using PlacedShape_t    = PlacedScaledShape;
  using UnplacedStruct_t = ScaledShapeStruct<Precision>;
  using UnplacedVolume_t = UnplacedScaledShape;

  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void Contains(UnplacedStruct_t const &unplaced,
                                                                    Vector3D<Real_v> const &point, bool &inside);

  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void Inside(UnplacedStruct_t const &unplaced,
                                                                  Vector3D<Real_v> const &point, Inside_t &inside);

  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void DistanceToIn(UnplacedStruct_t const &unplaced,
                                                                        Vector3D<Real_v> const &point,
                                                                        Vector3D<Real_v> const &direction,
                                                                        Real_v const &stepMax, Real_v &distance);

  /// @brief Compute the exit distance and optionally a supporting exit normal.
  /// @param unplaced Constituent and invertible diagonal scale.
  /// @param point Start point in the scaled local frame.
  /// @param direction Unit direction in the scaled local frame.
  /// @param stepMax Step limit forwarded in unscaled units, preserving the infinity sentinel.
  /// @param[out] distance Constituent exit distance converted to scaled units.
  /// @param[out] hitInfo Optional unit normal and supporting-plane flag.
  /// @details Supporting planes are preserved by an invertible scale. A false
  /// constituent flag is retained without writing the caller's normal storage.
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void DistanceToOut(UnplacedStruct_t const &unplaced,
                                                                         Vector3D<Real_v> const &point,
                                                                         Vector3D<Real_v> const &direction,
                                                                         Real_v const &stepMax, Real_v &distance,
                                                                         SurfaceHitView<Real_v> *hitInfo = nullptr);

  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SafetyToIn(UnplacedStruct_t const &unplaced,
                                                                      Vector3D<Real_v> const &point, Real_v &safety);

  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SafetyToOut(UnplacedStruct_t const &unplaced,
                                                                       Vector3D<Real_v> const &point, Real_v &safety);

  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void NormalKernel(UnplacedStruct_t const &unplaced,
                                                                        Vector3D<Real_v> const &point,
                                                                        Vector3D<Real_v> &normal, bool &valid);

}; // End struct ScaledShapeImplementation

// Implementations start here

template <typename Real_v>
VECCORE_ATT_HOST_DEVICE void ScaledShapeImplementation::Contains(UnplacedStruct_t const &unplaced,
                                                                 Vector3D<Real_v> const &point, bool &inside)
{

  // Transform local point to unscaled shape frame
  Vector3D<Real_v> ulocalPoint;
  unplaced.fScale.Transform(point, ulocalPoint);

  // Now call Contains for the unscaled shape
  inside = unplaced.fPlaced->Contains(ulocalPoint);
}

template <typename Real_v>
VECCORE_ATT_HOST_DEVICE void ScaledShapeImplementation::Inside(UnplacedStruct_t const &unplaced,
                                                               Vector3D<Real_v> const &point, Inside_t &inside)
{

  // Transform local point to unscaled shape frame
  Vector3D<Real_v> ulocalPoint;
  unplaced.fScale.Transform(point, ulocalPoint);

  // Now call Inside for the unscaled shape
  inside = unplaced.fPlaced->Inside(ulocalPoint);
}

template <typename Real_v>
VECCORE_ATT_HOST_DEVICE void ScaledShapeImplementation::DistanceToIn(UnplacedStruct_t const &unplaced,
                                                                     Vector3D<Real_v> const &point,
                                                                     Vector3D<Real_v> const &direction,
                                                                     Real_v const &stepMax, Real_v &distance)
{

  // Transform point, direction and stepMax to unscaled shape frame
  Vector3D<Real_v> ulocalPoint;
  unplaced.fScale.Transform(point, ulocalPoint);

  // Direction is un-normalized after scale transformation
  Vector3D<Real_v> ulocalDir;
  unplaced.fScale.Transform(direction, ulocalDir);
  ulocalDir.Normalize();

  // The maximum-distance sentinel is not a finite length to rescale. Saturate
  // large finite limits as well, so constituents never receive an overflowed limit.
  Real_v ustepMax = InfinityLength<Real_v>();
  if (stepMax < ustepMax) {
    const auto stepScale = unplaced.fScale.TransformDistance(Real_v(1.), direction);
    if (stepScale <= Real_v(1.) || stepMax < ustepMax / stepScale) ustepMax = stepMax * stepScale;
  }

  // Compute distance in unscaled system
  distance = unplaced.fPlaced->DistanceToIn(ulocalPoint, ulocalDir, ustepMax);

  // Convert finite distances back to master, preserving infinity sentinels.
  if (distance < InfinityLength<Real_v>()) distance = unplaced.fScale.InverseTransformDistance(distance, ulocalDir);
}

template <typename Real_v>
VECCORE_ATT_HOST_DEVICE void ScaledShapeImplementation::DistanceToOut(UnplacedStruct_t const &unplaced,
                                                                      Vector3D<Real_v> const &point,
                                                                      Vector3D<Real_v> const &direction,
                                                                      Real_v const &stepMax, Real_v &distance,
                                                                      SurfaceHitView<Real_v> *hitInfo)
{

  // Transform point, direction and stepMax to unscaled shape frame
  Vector3D<Real_v> ulocalPoint;
  unplaced.fScale.Transform(point, ulocalPoint);

  // Direction is un-normalized after scale transformation
  Vector3D<Real_v> ulocalDir;
  unplaced.fScale.Transform(direction, ulocalDir);
  ulocalDir.Normalize();

  // Preserve the unbounded-step sentinel and avoid overflowing finite limits.
  Real_v ustepMax = InfinityLength<Real_v>();
  if (stepMax < ustepMax) {
    const auto stepScale = unplaced.fScale.TransformDistance(Real_v(1.), direction);
    if (stepScale <= Real_v(1.) || stepMax < ustepMax / stepScale) ustepMax = stepMax * stepScale;
  }

  // An invertible scale preserves supporting planes; flag-only queries need
  // no normal transformation or reconstruction at the propagated hit point.
  if (!hitInfo) {
    distance = unplaced.fPlaced->DistanceToOut(ulocalPoint, ulocalDir, ustepMax);
  } else if (!hitInfo->WantsNormal()) {
    distance = unplaced.fPlaced->DistanceToOut(ulocalPoint, ulocalDir, ustepMax, hitInfo);
  } else {
    Vector3D<Real_v> ulocalNormal;
    SurfaceHitView<Real_v> ulocalHit{&ulocalNormal};
    distance = unplaced.fPlaced->DistanceToOut(ulocalPoint, ulocalDir, ustepMax, &ulocalHit);
    hitInfo->SetLeftBehind(ulocalHit.fLeftBehind);
    if (ulocalHit.fLeftBehind) hitInfo->SetNormal(unplaced.fScale.InverseTransformNormal(ulocalNormal));
  }

  // Convert finite distances back to master, preserving infinity sentinels.
  if (distance < InfinityLength<Real_v>()) distance = unplaced.fScale.InverseTransformDistance(distance, ulocalDir);
}

template <typename Real_v>
VECCORE_ATT_HOST_DEVICE void ScaledShapeImplementation::SafetyToIn(UnplacedStruct_t const &unplaced,
                                                                   Vector3D<Real_v> const &point, Real_v &safety)
{

  // Transform point to unscaled shape frame
  Vector3D<Real_v> ulocalPoint;
  unplaced.fScale.Transform(point, ulocalPoint);

  // Compute unscaled safety, then scale it.
  safety = unplaced.fPlaced->SafetyToIn(ulocalPoint);
  safety = unplaced.fScale.InverseTransformSafety(safety);
}

template <typename Real_v>
VECCORE_ATT_HOST_DEVICE void ScaledShapeImplementation::SafetyToOut(UnplacedStruct_t const &unplaced,
                                                                    Vector3D<Real_v> const &point, Real_v &safety)
{

  // Transform point to unscaled shape frame
  Vector3D<Real_v> ulocalPoint;
  unplaced.fScale.Transform(point, ulocalPoint);

  // Compute unscaled safety, then scale it.
  safety = unplaced.fPlaced->SafetyToOut(ulocalPoint);
  safety = unplaced.fScale.InverseTransformSafety(safety);
}

template <typename Real_v>
VECCORE_ATT_HOST_DEVICE void ScaledShapeImplementation::NormalKernel(UnplacedStruct_t const &unplaced,
                                                                     Vector3D<Real_v> const &point,
                                                                     Vector3D<Real_v> &normal, bool &valid)
{

  // Transform point to unscaled shape frame
  Vector3D<Real_v> ulocalPoint;
  unplaced.fScale.Transform(point, ulocalPoint);

  // Compute normal in unscaled frame
  Vector3D<Real_v> ulocalNorm;
  valid = unplaced.fPlaced->Normal(ulocalPoint, ulocalNorm);

  // Convert normal to scaled frame
  unplaced.fScale.InverseTransformNormal(ulocalNorm, normal);
}
} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif // VECGEOM_VOLUMES_KERNEL_SCALEDSHAPEIMPLEMENTATION_H_
