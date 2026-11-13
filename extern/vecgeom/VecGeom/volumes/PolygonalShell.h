/// @file PolygonalShell.h
/// @brief Z-aligned side-shell helper for simple extruded polygons.

#ifndef VECGEOM_POLYGONAL_SHELL_H
#define VECGEOM_POLYGONAL_SHELL_H

#include "VecGeom/base/Global.h"
#include "VecGeom/volumes/PlanarPolygon.h"

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE(class PolygonalShell;);
VECGEOM_DEVICE_DECLARE_CONV(class, PolygonalShell);

inline namespace VECGEOM_IMPL_NAMESPACE {

/// @brief Side-shell representation for a polygon extruded between two z planes.
/// @details Each polygon edge defines one z-aligned rectangular side surface.
/// The shell provides side-only distance helpers used by `SExtruImplementation`;
/// cap handling stays in the SExtru kernel.
class PolygonalShell : AlignedBase {

private:
  // the polygon (with friend access)
  PlanarPolygon fPolygon;
  Precision fLowerZ;         // lower z plane
  Precision fUpperZ;         // upper z plane
  Vector<int> fSideSupports; ///< Whether each outward side plane supports the complete polygon

  /// @brief Clip an outward side-plane hit to its segment and forward material ownership.
  /// @details The caller establishes plane membership. Unit side tangents test
  /// segment bounds in length units, without repeating a cancellation-prone
  /// plane test after propagation. At a concave vertex the inward half-spaces form a union.
  /// Leaving this side is not an exit when the ray enters the adjoining side.
  /// Convex vertices and the interior of a side retain their ordinary ownership.
  /// @param side Side index with an outward-directed crossing.
  /// @param x,y Hit coordinates in the side's plane tolerance band.
  /// @param dir Local ray direction.
  /// @return Whether this finite side owns the exit.
  template <typename Real_v>
  VECCORE_ATT_HOST_DEVICE bool IsExitOnSide(size_t side, Real_v x, Real_v y, Vector3D<Real_v> const &dir) const
  {
    const size_t count = fPolygon.GetNVertices();
    const Real_v along =
        (x - fPolygon.fVertices.x()[side]) * fPolygon.fB[side] - (y - fPolygon.fVertices.y()[side]) * fPolygon.fA[side];
    const Real_v beyond =
        (x - fPolygon.fShiftedXJ[side]) * fPolygon.fB[side] - (y - fPolygon.fShiftedYJ[side]) * fPolygon.fA[side];
    if (along < -kToleranceDist<Real_v> || beyond > kToleranceDist<Real_v>) return false;
    for (int end = 0; end < 2; ++end) {
      if (Abs(end ? along : beyond) >= kToleranceDist<Real_v>) continue;
      const size_t vertex = end ? side : (side ? side - 1 : count - 1);
      const size_t next   = vertex + 1 == count ? 0 : vertex + 1;
      if (fPolygon.fA[vertex] * fPolygon.fB[next] - fPolygon.fB[vertex] * fPolygon.fA[next] <= 0.) continue;
      const size_t other = end ? next : vertex;
      if (fPolygon.fA[other] * dir.x() + fPolygon.fB[other] * dir.y() > kToleranceDist<Real_v>) return false;
    }
    return true;
  }

  /// @brief Cache whether each side plane leaves the complete polygon behind.
  VECCORE_ATT_HOST_DEVICE
  void CacheSideSupport()
  {
    const size_t sideCount = fPolygon.GetNVertices();
    fSideSupports.resize(sideCount, 1);
    if (fPolygon.IsConvex()) return;

    auto const &vertices = fPolygon.GetVertices();
    for (size_t side = 0; side < sideCount; ++side) {
      bool supports = true;
      for (size_t vertex = 0; vertex < sideCount; ++vertex) {
        const Precision distance =
            -(fPolygon.fA[side] * vertices.x()[vertex] + fPolygon.fB[side] * vertices.y()[vertex] + fPolygon.fD[side]);
        if (distance > kToleranceDist<Precision>) {
          supports = false;
          break;
        }
      }
      fSideSupports[side] = static_cast<int>(supports);
    }
  }

  friend class SimpleExtruPolygon;
  friend struct SExtruImplementation;
  friend class UnplacedSExtruVolume;

public:
  VECCORE_ATT_HOST_DEVICE
  PolygonalShell() : fPolygon() {}

  VECCORE_ATT_HOST_DEVICE
  PolygonalShell(int nvertices, Precision *x, Precision *y, Precision lowerz, Precision upperz)
      : fPolygon(nvertices, x, y), fLowerZ(lowerz), fUpperZ(upperz)
  {
    CacheSideSupport();
  }

  VECCORE_ATT_HOST_DEVICE
  void Init(int nvertices, Precision *x, Precision *y, Precision lowerz, Precision upperz)
  {
    fPolygon.Init(nvertices, x, y);
    fLowerZ = lowerz;
    fUpperZ = upperz;
    CacheSideSupport();
  }

  // the area of the shell ( does not include the area of the planar polygon )
  VECCORE_ATT_HOST_DEVICE
  Precision SurfaceArea() const
  {
    const auto kS = fPolygon.fVertices.size();
    Precision area(0.);
    for (size_t i = 0; i < kS; ++i) {
      area += Sqrt(fPolygon.fLengthSqr[i]);
    }
    return area * (fUpperZ - fLowerZ);
  }

  VECCORE_ATT_HOST_DEVICE
  PlanarPolygon const &GetPolygon() const { return fPolygon; }

  VECCORE_ATT_HOST_DEVICE
  Precision GetLowerZ() const { return fLowerZ; }

  VECCORE_ATT_HOST_DEVICE
  Precision GetUpperZ() const { return fUpperZ; }

  /// @brief Return whether a side's outward plane supports the complete polygon.
  /// @param side Polygon edge and side-surface index.
  /// @return True if all polygon vertices are behind its outward plane.
  VECCORE_ATT_HOST_DEVICE
  bool SideSupports(size_t side) const { return fSideSupports[side] != 0; }

  /// @brief Return the cached outward unit normal of a side surface.
  /// @param side Polygon edge and side-surface index.
  /// @return Outward unit normal, with a zero z component.
  template <typename Real_v>
  VECCORE_ATT_HOST_DEVICE Vector3D<Real_v> GetSideNormal(size_t side) const
  {
    return Vector3D<Real_v>(-fPolygon.fA[side], -fPolygon.fB[side], Real_v(0.));
  }

  template <typename Real_v>
  VECCORE_ATT_HOST_DEVICE void Extent(Vector3D<Real_v> &aMin, Vector3D<Real_v> &aMax) const
  {
    aMin[0] = Real_v(fPolygon.GetMinX());
    aMin[1] = Real_v(fPolygon.GetMinY());
    aMin[2] = Real_v(fLowerZ);

    aMax[0] = Real_v(fPolygon.GetMaxX());
    aMax[1] = Real_v(fPolygon.GetMaxY());
    aMax[2] = Real_v(fUpperZ);
  }

  template <typename Real_v>
  VECCORE_ATT_HOST_DEVICE Real_v DistanceToIn(Vector3D<Real_v> const &point, Vector3D<Real_v> const &dir) const
  {
    return fPolygon.IsConvex() ? DistanceToInConvex(point, dir) : DistanceToInConcave(point, dir);
  }

  /// @brief Compute the first side-shell entry for a convex polygon.
  /// @details Intersects candidate side planes, then filters by z range and
  /// segment ownership. Convex shells can return as soon as the first valid
  /// side hit is found.
  /// @param point Local start point.
  /// @param dir Unit local direction.
  /// @return Side entry distance or `kInfLength` when no side is hit.
  template <typename Real_v>
  VECCORE_ATT_HOST_DEVICE Real_v DistanceToInConvex(Vector3D<Real_v> const &point, Vector3D<Real_v> const &dir) const
  {
    bool done = false;
    Real_v result(kInfLength);
    const auto S = fPolygon.fVertices.size();

    for (size_t i = 0; i < S; ++i) { // side/rectangle index
      // approaching from right side?
      // under the assumption that surface normals points "inwards"
      const Real_v proj      = fPolygon.fA[i] * dir.x() + fPolygon.fB[i] * dir.y();
      const bool sidecorrect = proj >= -kTolerance;
      if (!sidecorrect) continue;

      // the distance to the plane (specialized for fNormalsZ == 0)
      const Real_v pdist = fPolygon.fA[i] * point.x() + fPolygon.fB[i] * point.y() + fPolygon.fD[i];

      const bool moving_away = pdist > kTolerance;
      if (moving_away) continue;

      const Real_v dist = -pdist / NonZero(proj);

      // propagate to plane (first just z)
      const Real_v zInters(point.z() + dist * dir.z());
      const Real_v tolerance = kToleranceDist<Real_v>;
      const bool zRangeOk    = (zInters <= fUpperZ + tolerance) && (zInters >= fLowerZ - tolerance);
      if (zRangeOk) {
        // check intersection with rest of rectangle
        const Real_v xInters(point.x() + dist * dir.x());
        const Real_v yInters(point.y() + dist * dir.y());

        // we could already check if intersection within the known extent
        const bool intersects = fPolygon.OnSegment<Real_v, Precision, bool>(i, xInters, yInters);

        if (!done && intersects) {
          result = Max(dist, Real_v(0.));
          done   = true;
        }
      }
      if (done) return result;
    }
    return result;
  }

  /// @brief Compute the nearest side-shell entry for a concave polygon.
  /// @details Concave shells may expose multiple side planes to the ray, so all
  /// accepted plane/segment candidates are reduced to the nearest non-negative
  /// entry distance.
  /// @param point Local start point.
  /// @param dir Unit local direction.
  /// @return Side entry distance or `kInfLength` when no side is hit.
  template <typename Real_v>
  VECCORE_ATT_HOST_DEVICE Real_v DistanceToInConcave(Vector3D<Real_v> const &point, Vector3D<Real_v> const &dir) const
  {
    Real_v result(kInfLength);
    const auto S = fPolygon.fVertices.size();

    for (size_t i = 0; i < S; ++i) { // side/rectangle index
      // approaching from right side?
      // under the assumption that surface normals points "inwards"
      const Real_v proj = fPolygon.fA[i] * dir.x() + fPolygon.fB[i] * dir.y();
      // A negative inward projection is an exit, even at shallow incidence.
      const bool sidecorrect = proj >= Real_v(0.);
      if (!sidecorrect) continue;

      // the distance to the plane (specialized for fNormalsZ == 0)
      const Real_v pdist = fPolygon.fA[i] * point.x() + fPolygon.fB[i] * point.y() + fPolygon.fD[i];

      const bool moving_away = pdist > kTolerance;
      if (moving_away) continue;

      const Real_v dist = -pdist / NonZero(proj);

      // propagate to plane (first just z)
      const Real_v zInters(point.z() + dist * dir.z());
      const Real_v tolerance = kToleranceDist<Real_v>;
      const bool zRangeOk    = (zInters <= fUpperZ + tolerance) && (zInters >= fLowerZ - tolerance);
      if (zRangeOk) {
        // check intersection with rest of rectangle
        const Real_v xInters(point.x() + dist * dir.x());
        const Real_v yInters(point.y() + dist * dir.y());

        // we could already check if intersection within the known extent
        const bool intersects = fPolygon.OnSegment<Real_v, Precision, bool>(i, xInters, yInters);

        if (intersects) result = Min(Max(dist, Real_v(0.)), result);
      }
    }
    return result;
  }

  // -- DistanceToOut --

  template <typename Real_v>
  VECCORE_ATT_HOST_DEVICE Real_v DistanceToOut(Vector3D<Real_v> const &point, Vector3D<Real_v> const &dir,
                                               int *hitSide = nullptr) const
  {
    return fPolygon.IsConvex() ? DistanceToOutConvex(point, dir, hitSide) : DistanceToOutConcave(point, dir, hitSide);
  }

  /// @brief Compute the first side-shell exit for a convex polygon.
  /// @details This is the exit counterpart of `DistanceToInConvex`: candidate
  /// side-plane intersections are filtered by z range and segment ownership.
  /// Tolerated side-surface starts moving outward return zero.
  /// @param point Local start point.
  /// @param dir Unit local direction.
  /// @param[out] hitSide Optional selected side-surface index.
  /// @return Side exit distance or `kInfLength` when no side is hit.
  template <typename Real_v>
  VECCORE_ATT_HOST_DEVICE Real_v DistanceToOutConvex(Vector3D<Real_v> const &point, Vector3D<Real_v> const &dir,
                                                     int *hitSide = nullptr) const
  {
    if (hitSide) *hitSide = -1;
    bool done = false;
    Real_v result(kInfLength);
    const auto S = fPolygon.fVertices.size();

    for (size_t i = 0; i < S; ++i) { // side/rectangle index
      // approaching from right side?
      // under the assumption that surface normals points "inwards"
      const Real_v proj = fPolygon.fA[i] * dir.x() + fPolygon.fB[i] * dir.y();

      // the distance to the plane (specialized for fNormalsZ == 0)
      const Real_v pdist = fPolygon.fA[i] * point.x() + fPolygon.fB[i] * point.y() + fPolygon.fD[i];

      // Tolerated side-surface starts moving out are immediate exits.
      const bool surfaceExit = (Abs(pdist) <= kToleranceDist<Real_v>) && (proj < Real_v(0.)) &&
                               (point.z() <= fUpperZ + kToleranceDist<Real_v>) &&
                               (point.z() >= fLowerZ - kToleranceDist<Real_v>) &&
                               fPolygon.OnSegment<Real_v, Precision, bool>(i, point.x(), point.y());
      if (!done && surfaceExit) {
        result = Real_v(0.);
        done   = true;
        if (hitSide) *hitSide = static_cast<int>(i);
      }
      if (done) return result;

      const bool sidecorrect = proj <= -kTolerance;
      if (!sidecorrect) continue;

      const bool moving_away = pdist < -kTolerance;
      if (moving_away) continue;

      const Real_v dist = -pdist / NonZero(proj);

      // propagate to plane (first just z)
      const Real_v zInters(point.z() + dist * dir.z());
      const Real_v tolerance = kToleranceDist<Real_v>;
      const bool zRangeOk =
          (zInters <= fUpperZ + tolerance) && (zInters >= fLowerZ - tolerance) && sidecorrect && !moving_away;
      if (zRangeOk) {
        // check intersection with rest of rectangle
        const Real_v xInters(point.x() + dist * dir.x());
        const Real_v yInters(point.y() + dist * dir.y());

        // we could already check if intersection within the known extent
        const bool intersects =
            fPolygon.OnSegment<Real_v, Precision, bool>(i, xInters, yInters) && (dist >= -Real_v(kTolerance));

        if (!done && intersects) {
          result = dist;
          done   = true;
          if (hitSide) *hitSide = static_cast<int>(i);
        }
      }
      if (done) return result;
    }
    return result;
  }

  /// @brief Compute the nearest side-shell exit for a concave polygon.
  /// @details Concave shells may have multiple valid side candidates, so the
  /// nearest accepted plane/segment hit is selected. Tolerated side-surface
  /// starts moving outward return zero.
  /// @param point Local start point.
  /// @param dir Unit local direction.
  /// @param[out] hitSide Optional selected side-surface index.
  /// @return Side exit distance or `kInfLength` when no side is hit.
  template <typename Real_v>
  VECCORE_ATT_HOST_DEVICE Real_v DistanceToOutConcave(Vector3D<Real_v> const &point, Vector3D<Real_v> const &dir,
                                                      int *hitSide = nullptr) const
  {
    if (hitSide) *hitSide = -1;
    Real_v result(kInfLength);
    const auto S = fPolygon.fVertices.size();

    for (size_t i = 0; i < S; ++i) { // side/rectangle index
      // approaching from right side?
      // under the assumption that surface normals points "inwards"
      const Real_v proj = fPolygon.fA[i] * dir.x() + fPolygon.fB[i] * dir.y();

      // the distance to the plane (specialized for fNormalsZ == 0)
      const Real_v pdist = fPolygon.fA[i] * point.x() + fPolygon.fB[i] * point.y() + fPolygon.fD[i];

      // Tolerated side-surface starts moving out are immediate exits.
      const bool surfaceExit = (Abs(pdist) <= kToleranceDist<Real_v>) && (proj < Real_v(0.)) &&
                               (point.z() <= fUpperZ + kToleranceDist<Real_v>) &&
                               (point.z() >= fLowerZ - kToleranceDist<Real_v>) &&
                               IsExitOnSide(i, point.x(), point.y(), dir);
      if (surfaceExit) {
        if (hitSide) *hitSide = static_cast<int>(i);
        return Real_v(0.);
      }

      const bool sidecorrect = proj < -kTolerance;
      if (!sidecorrect) continue;

      const bool moving_away = pdist < -kTolerance;
      if (moving_away) continue;

      const Real_v dist = -pdist / NonZero(proj);

      // propagate to plane (first just z)
      const Real_v zInters(point.z() + dist * dir.z());
      const Real_v tolerance = kToleranceDist<Real_v>;
      const bool zRangeOk =
          (zInters <= fUpperZ + tolerance) && (zInters >= fLowerZ - tolerance) && sidecorrect && !moving_away;
      if (zRangeOk) {
        // check intersection with rest of rectangle
        const Real_v xInters(point.x() + dist * dir.x());
        const Real_v yInters(point.y() + dist * dir.y());

        // we could already check if intersection within the known extent
        const bool intersects = (dist >= -Real_v(kTolerance)) && IsExitOnSide(i, xInters, yInters, dir);

        // Accepted negative roots are tolerated surface starts, not backward exits.
        if (intersects) {
          const Real_v candidate = Max(dist, Real_v(0.));
          if (hitSide && candidate < result) *hitSide = static_cast<int>(i);
          result = Min(candidate, result);
        }
      }
    }
    return result;
  }

}; // end class

#define SPECIALIZATION
#ifdef SPECIALIZATION
// template specialization for Distance functions
template <>
VECCORE_ATT_HOST_DEVICE inline Precision PolygonalShell::DistanceToOutConvex(Vector3D<Precision> const &point,
                                                                             Vector3D<Precision> const &dir,
                                                                             int *hitSide) const
{
  if (hitSide) *hitSide = -1;
  Precision dz         = 0.5 * (fUpperZ - fLowerZ);
  Precision pz         = point.z() - 0.5 * (fLowerZ + fUpperZ);
  const Precision safz = Abs(pz) - dz;
  if (safz > kTolerance) return -kTolerance;

  Precision vz        = dir.z();
  Precision tmax      = kInfLength;
  int selectedSurface = -1;
  Precision minOffset = 0.; // Reduce inward offsets; never publish an outside hit.
  const auto S        = fPolygon.fVertices.size();
  if (Abs(vz) > kTolerance) {
    tmax            = (CopySign(dz, vz) - pz) / vz;
    selectedSurface = static_cast<int>(S) + (vz > 0. ? 1 : 0);
  }
  for (size_t i = 0; i < S; ++i) { // side/rectangle index

    // Stored planes point inward: negative offsets are outside, negative projections leave.
    const Precision proj  = fPolygon.fA[i] * dir.x() + fPolygon.fB[i] * dir.y();
    const Precision pdist = fPolygon.fA[i] * point.x() + fPolygon.fB[i] * point.y() + fPolygon.fD[i];
    minOffset             = Min(pdist, minOffset);
    if (proj < -kTolerance) {
      const Precision dist = -pdist / proj;
      if (tmax > dist) {
        tmax            = dist;
        selectedSurface = static_cast<int>(i);
      }
    }
  }
  if (minOffset < -kTolerance) return -kTolerance;
  if (hitSide) *hitSide = selectedSurface;
  // Accepted surface exits may be slightly negative from plane arithmetic; the
  // navigation convention is an immediate zero step, not a backwards step.
  return Max(tmax, Precision(0.));
}

// template specialization for Distance functions
template <>
VECCORE_ATT_HOST_DEVICE inline Precision PolygonalShell::DistanceToInConvex(Vector3D<Precision> const &point,
                                                                            Vector3D<Precision> const &dir) const
{
  Precision dz = 0.5 * (fUpperZ - fLowerZ);
  Precision pz = point.z() - 0.5 * (fLowerZ + fUpperZ);
  if ((Abs(pz) - dz) > -kTolerance && pz * dir.z() >= 0) return kInfLength;
  const Precision invz = -1. / NonZero(dir.z());
  const Precision ddz  = (invz < 0) ? dz : -dz;
  Precision tmin       = (pz + ddz) * invz;
  Precision tmax       = (pz - ddz) * invz;
  const auto S         = fPolygon.fVertices.size();
  for (size_t i = 0; i < S; ++i) { // side/rectangle index

    const Precision proj = -(fPolygon.fA[i] * dir.x() + fPolygon.fB[i] * dir.y());
    // the distance to the plane (specialized for fNormalsZ == 0)
    const Precision pdist   = -(fPolygon.fA[i] * point.x() + fPolygon.fB[i] * point.y() + fPolygon.fD[i]);
    const bool side_correct = pdist > -kTolerance;
    // A parallel plane rejects an outside start but does not constrain the interval from an inside start.
    if (proj == 0.) {
      if (side_correct) return kInfLength;
      continue;
    }
    const bool moving_away = proj > 0;
    if (side_correct) {
      if (moving_away) return kInfLength;
      const Precision dist = -pdist / NonZero(proj);
      if (dist > tmin) tmin = dist;
    } else if (moving_away) {
      const Precision dist = -pdist / NonZero(proj);
      if (dist < tmax) tmax = dist;
    }
  }
  if (tmax < tmin + kTolerance) return kInfLength;
  return tmin;
}

#endif

} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif
