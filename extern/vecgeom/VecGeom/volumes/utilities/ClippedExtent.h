// This file is part of VecGeom and is distributed under the
// conditions in the file LICENSE.txt in the top directory.
// For the full list of authors see CONTRIBUTORS.txt and `git log`.

/// \brief Utilities for host-side clipped extent queries.
/// \file volumes/utilities/ClippedExtent.h

#ifndef VECGEOM_VOLUMES_UTILITIES_CLIPPEDEXTENT_H_
#define VECGEOM_VOLUMES_UTILITIES_CLIPPEDEXTENT_H_

#include "VecGeom/base/Global.h"
#include "VecGeom/base/Transformation3D.h"
#include "VecGeom/base/Vector3D.h"

namespace vecgeom {

/// @brief Coordinate axis selector for clipped extent queries.
enum class Axis { kX = 0, kY = 1, kZ = 2 };

/// @brief Convert an axis selector to a `Vector3D` component index.
/// @param axis Axis selector.
/// @return Component index in `[0, 2]`.
inline int AxisIndex(Axis axis) { return static_cast<int>(axis); }

/// @brief Convert a component index to an axis selector.
/// @param axis Component index in `[0, 2]`.
/// @return Axis selector with the same numeric value.
inline Axis AxisFromIndex(int axis) { return static_cast<Axis>(axis); }

/// @brief Return the lower-case coordinate name of an axis.
/// @param axis Axis selector.
/// @return `"x"`, `"y"`, or `"z"`.
inline const char *AxisName(Axis axis) { return axis == Axis::kX ? "x" : (axis == Axis::kY ? "y" : "z"); }

/// @brief One-dimensional conservative extent interval.
struct ExtentInterval {
  /// True when `min` and `max` describe a non-empty interval.
  bool exists = false;
  /// Lower interval bound when `exists` is true.
  Precision min = kInfLength;
  /// Upper interval bound when `exists` is true.
  Precision max = -kInfLength;

  ExtentInterval() = default;

  /// @brief Construct an interval from explicit bounds.
  /// @param min_value Lower interval bound.
  /// @param max_value Upper interval bound.
  ExtentInterval(Precision min_value, Precision max_value)
      : exists(min_value <= max_value), min(min_value), max(max_value)
  {
  }
};

/// @brief Conservative axis-aligned box returned by clipped extent queries.
struct ClippedBox {
  /// True when `min` and `max` describe a non-empty box.
  bool exists = false;
  /// Lower box corner when `exists` is true.
  Vector3D<Precision> min{kInfLength, kInfLength, kInfLength};
  /// Upper box corner when `exists` is true.
  Vector3D<Precision> max{-kInfLength, -kInfLength, -kInfLength};

  ClippedBox() = default;

  /// @brief Construct a box from explicit lower and upper corners.
  /// @param min_value Lower box corner.
  /// @param max_value Upper box corner.
  ClippedBox(Vector3D<Precision> const &min_value, Vector3D<Precision> const &max_value)
      : exists(min_value.x() <= max_value.x() && min_value.y() <= max_value.y() && min_value.z() <= max_value.z()),
        min(min_value), max(max_value)
  {
  }

  /// @brief Return this box projected onto one axis.
  /// @param axis Projection axis.
  /// @return Existing interval when this box exists, otherwise an empty interval.
  ExtentInterval Interval(Axis axis) const
  {
    const int index = AxisIndex(axis);
    return exists ? ExtentInterval(min[index], max[index]) : ExtentInterval();
  }

  /// @brief Expand this box to include a point.
  /// @param point Point to include.
  void Include(Vector3D<Precision> const &point)
  {
    if (!exists) {
      min    = point;
      max    = point;
      exists = true;
      return;
    }
    for (int axis = 0; axis < 3; ++axis) {
      if (point[axis] < min[axis]) min[axis] = point[axis];
      if (point[axis] > max[axis]) max[axis] = point[axis];
    }
  }

  /// @brief Expand this box to include another box.
  /// @param other Box to include; empty boxes are ignored.
  void Include(ClippedBox const &other)
  {
    if (!other.exists) return;
    Include(other.min);
    Include(other.max);
  }

  /// @brief Enlarge all existing box bounds by a symmetric margin.
  /// @param margin Non-negative margin to subtract from `min` and add to `max`.
  void Inflate(Precision margin)
  {
    if (!exists || margin <= Precision(0.)) return;
    min -= Vector3D<Precision>(margin);
    max += Vector3D<Precision>(margin);
  }

  /// @brief Test whether a point is inside the box within tolerance.
  /// @param point Point in the same frame as the box.
  /// @param tolerance Symmetric tolerance applied to all faces.
  /// @return True when the box exists and contains the point within tolerance.
  bool Contains(Vector3D<Precision> const &point, Precision tolerance = Precision(0.)) const
  {
    if (!exists) return false;
    return point.x() >= min.x() - tolerance && point.x() <= max.x() + tolerance && point.y() >= min.y() - tolerance &&
           point.y() <= max.y() + tolerance && point.z() >= min.z() - tolerance && point.z() <= max.z() + tolerance;
  }
};

/// @brief Axis-aligned mother-frame clipping limits for extent queries.
class VoxelLimits {
public:
  /// @brief Construct unlimited, non-empty clipping limits.
  VoxelLimits()
  {
    for (int axis = 0; axis < 3; ++axis) {
      fMin[axis]     = -kInfLength;
      fMax[axis]     = kInfLength;
      fLimited[axis] = false;
    }
  }

  /// @brief Intersect the current limits with one axis interval.
  /// @param axis Axis to restrict.
  /// @param min_value Lower accepted coordinate on `axis`.
  /// @param max_value Upper accepted coordinate on `axis`.
  void AddLimit(Axis axis, Precision min_value, Precision max_value)
  {
    const int index = AxisIndex(axis);
    if (min_value > max_value) {
      fEmpty = true;
      return;
    }
    if (fLimited[index]) {
      if (min_value > fMin[index]) fMin[index] = min_value;
      if (max_value < fMax[index]) fMax[index] = max_value;
    } else {
      fMin[index]     = min_value;
      fMax[index]     = max_value;
      fLimited[index] = true;
    }
    if (fMin[index] > fMax[index]) fEmpty = true;
  }

  /// @brief Query whether one axis has finite limits.
  /// @param axis Axis to query.
  /// @return True when `axis` has an active finite interval.
  bool IsLimited(Axis axis) const { return fLimited[AxisIndex(axis)]; }

  /// @brief Query whether any axis has finite limits.
  /// @return True when at least one axis is limited.
  bool IsLimited() const { return fLimited[0] || fLimited[1] || fLimited[2]; }

  /// @brief Query whether previous limit intersections made the voxel empty.
  /// @return True when no point can satisfy the limits.
  bool IsEmpty() const { return fEmpty; }

  /// @brief Lower coordinate limit on one axis.
  /// @param axis Axis to query.
  /// @return Lower limit, or `-kInfLength` for an unlimited axis.
  Precision Min(Axis axis) const { return fMin[AxisIndex(axis)]; }

  /// @brief Upper coordinate limit on one axis.
  /// @param axis Axis to query.
  /// @return Upper limit, or `kInfLength` for an unlimited axis.
  Precision Max(Axis axis) const { return fMax[AxisIndex(axis)]; }

  /// @brief Test whether a point satisfies all active limits.
  /// @param point Mother-frame point.
  /// @param tolerance Symmetric tolerance applied to active limits.
  /// @return True when the limits are non-empty and the point is inside.
  bool Inside(Vector3D<Precision> const &point, Precision tolerance = Precision(0.)) const
  {
    if (fEmpty) return false;
    for (int axis = 0; axis < 3; ++axis) {
      if (!fLimited[axis]) continue;
      if (point[axis] < fMin[axis] - tolerance || point[axis] > fMax[axis] + tolerance) return false;
    }
    return true;
  }

  /// @brief Clip a line segment to the active limits.
  /// @param[in,out] start Segment start; replaced by the clipped start on success.
  /// @param[in,out] end Segment end; replaced by the clipped end on success.
  /// @return True when the segment intersects the limits.
  bool ClipSegment(Vector3D<Precision> &start, Vector3D<Precision> &end) const
  {
    if (fEmpty) return false;
    const Vector3D<Precision> original_start = start;
    const Vector3D<Precision> delta          = end - start;
    Precision tmin                           = 0.;
    Precision tmax                           = 1.;

    for (int axis = 0; axis < 3; ++axis) {
      if (!fLimited[axis]) continue;
      const Precision direction = delta[axis];
      if (Abs(direction) <= kTolerance) {
        if (start[axis] < fMin[axis] || start[axis] > fMax[axis]) return false;
        continue;
      }
      Precision lo = (fMin[axis] - start[axis]) / direction;
      Precision hi = (fMax[axis] - start[axis]) / direction;
      if (lo > hi) {
        const Precision tmp = lo;
        lo                  = hi;
        hi                  = tmp;
      }
      if (lo > tmin) tmin = lo;
      if (hi < tmax) tmax = hi;
      if (tmin > tmax) return false;
    }

    start = original_start + tmin * delta;
    end   = original_start + tmax * delta;
    return true;
  }

private:
  Vector3D<Precision> fMin;
  Vector3D<Precision> fMax;
  bool fLimited[3] = {false, false, false};
  bool fEmpty      = false;
};

/// @brief Intersect a conservative box with voxel limits.
/// @param[in,out] box Box to clip; marked empty if there is no overlap.
/// @param limits Mother-frame limits.
/// @param margin Symmetric margin applied to limits before clipping.
/// @return True when the clipped box remains non-empty.
inline bool IntersectBoxWithLimits(ClippedBox &box, VoxelLimits const &limits, Precision margin = Precision(0.))
{
  if (!box.exists || limits.IsEmpty()) {
    box.exists = false;
    return false;
  }
  for (int axis = 0; axis < 3; ++axis) {
    const Axis clipped_axis = AxisFromIndex(axis);
    if (!limits.IsLimited(clipped_axis)) continue;
    const Precision lower = limits.Min(clipped_axis) - margin;
    const Precision upper = limits.Max(clipped_axis) + margin;
    if (box.max[axis] < lower || box.min[axis] > upper) {
      box.exists = false;
      return false;
    }
    if (box.min[axis] < lower) box.min[axis] = lower;
    if (box.max[axis] > upper) box.max[axis] = upper;
  }
  return box.exists;
}

/// @brief Transform a local AABB to a padded mother-frame AABB.
/// @param local_min Local lower corner.
/// @param local_max Local upper corner.
/// @param transform Local-to-mother transform used by the existing placement path.
/// @param padding Symmetric padding applied after transformation.
/// @return Conservative mother-frame AABB.
inline ClippedBox TransformLocalAABBToMother(Vector3D<Precision> const &local_min, Vector3D<Precision> const &local_max,
                                             Transformation3D const &transform, Precision padding = kAABBPadding)
{
  const Vector3D<Precision> delta = local_max - local_min;
  ClippedBox box;
  for (int x = 0; x <= 1; ++x) {
    for (int y = 0; y <= 1; ++y) {
      for (int z = 0; z <= 1; ++z) {
        Vector3D<Precision> corner(local_min.x() + x * delta.x(), local_min.y() + y * delta.y(),
                                   local_min.z() + z * delta.z());
        box.Include(transform.InverseTransform(corner));
      }
    }
  }
  box.Inflate(padding);
  return box;
}

/// @brief Transform a local AABB using a pointer-valued transform.
/// @param local_min Local lower corner.
/// @param local_max Local upper corner.
/// @param transform Local-to-mother transform pointer.
/// @param padding Symmetric padding applied after transformation.
/// @return Conservative mother-frame AABB.
inline ClippedBox TransformLocalAABBToMother(Vector3D<Precision> const &local_min, Vector3D<Precision> const &local_max,
                                             Transformation3D const *transform, Precision padding = kAABBPadding)
{
  return TransformLocalAABBToMother(local_min, local_max, *transform, padding);
}

} // namespace vecgeom

#endif // VECGEOM_VOLUMES_UTILITIES_CLIPPEDEXTENT_H_
