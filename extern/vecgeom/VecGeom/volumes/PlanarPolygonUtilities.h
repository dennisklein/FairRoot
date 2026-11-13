/// @file PlanarPolygonUtilities.h
/// @brief Geometry predicates shared by PlanarPolygon triangulation users.

#ifndef VECGEOM_PLANAR_POLYGON_UTILITIES_H
#define VECGEOM_PLANAR_POLYGON_UTILITIES_H

#include "VecGeom/base/Global.h"
#include "VecGeom/base/Vector3D.h"
#include "VecGeom/volumes/PlanarPolygon.h"

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {
namespace PlanarPolygonUtilities {

/// @brief Return a base polygon vertex as an XY-plane vector.
VECGEOM_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
Vector3D<Precision> Vertex2D(PlanarPolygon const &polygon, size_t index)
{
  auto const &vertices = polygon.GetVertices();
  return Vector3D<Precision>(vertices.x()[index], vertices.y()[index], 0.);
}

/// @brief Check whether a polygon vertex is on the line through two other vertices.
/// @details Uses a length-scaled squared cross product to avoid square roots.
VECGEOM_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
bool IsSameLine(PlanarPolygon const &polygon, size_t point, size_t first, size_t second)
{
  const auto *x        = polygon.GetVertices().x();
  const auto *y        = polygon.GetVertices().y();
  const Precision dx   = x[second] - x[first];
  const Precision dy   = y[second] - y[first];
  const Precision tol  = static_cast<Precision>(0.5) * kTolerance;
  const Precision tol2 = tol * tol;

  const Precision relx = x[point] - x[first];
  const Precision rely = y[point] - y[first];
  const Precision len2 = dx * dx + dy * dy;
  if (len2 <= tol2) return relx * relx + rely * rely <= tol2;

  const Precision cross = dx * rely - dy * relx;
  return cross * cross <= len2 * tol2;
}

/// @brief Check whether a polygon vertex is on a segment within tolerance.
VECGEOM_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
bool IsSameLineSegment(PlanarPolygon const &polygon, size_t point, size_t first, size_t second)
{
  const auto *x        = polygon.GetVertices().x();
  const auto *y        = polygon.GetVertices().y();
  const Precision tol  = static_cast<Precision>(0.5) * kTolerance;
  const Precision px   = x[point];
  const Precision py   = y[point];
  const Precision minx = Min(x[first], x[second]) - tol;
  const Precision maxx = Max(x[first], x[second]) + tol;
  const Precision miny = Min(y[first], y[second]) - tol;
  const Precision maxy = Max(y[first], y[second]) + tol;
  if (px < minx || px > maxx || py < miny || py > maxy) return false;
  return IsSameLine(polygon, point, first, second);
}

/// @brief Check whether two vertices are strictly on the same side of an edge.
VECGEOM_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
bool IsSameSide(PlanarPolygon const &polygon, size_t lhs, size_t rhs, size_t first, size_t second)
{
  const auto *x = polygon.GetVertices().x();
  const auto *y = polygon.GetVertices().y();

  return ((x[lhs] - x[first]) * (y[second] - y[first]) - (x[second] - x[first]) * (y[lhs] - y[first])) *
             ((x[rhs] - x[first]) * (y[second] - y[first]) - (x[second] - x[first]) * (y[rhs] - y[first])) >
         0.;
}

/// @brief Check whether a vertex lies inside or on a candidate triangle.
VECGEOM_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
bool IsPointInsideTriangle(PlanarPolygon const &polygon, size_t first, size_t second, size_t third, size_t point)
{
  const auto *x = polygon.GetVertices().x();
  const auto *y = polygon.GetVertices().y();

  if ((x[point] < x[first] && x[point] < x[second] && x[point] < x[third]) ||
      (x[point] > x[first] && x[point] > x[second] && x[point] > x[third]) ||
      (y[point] < y[first] && y[point] < y[second] && y[point] < y[third]) ||
      (y[point] > y[first] && y[point] > y[second] && y[point] > y[third]))
    return false;

  const bool inside  = IsSameSide(polygon, point, first, second, third) &&
                       IsSameSide(polygon, point, second, first, third) &&
                       IsSameSide(polygon, point, third, first, second);
  const bool on_edge = IsSameLineSegment(polygon, point, first, second) ||
                       IsSameLineSegment(polygon, point, second, third) ||
                       IsSameLineSegment(polygon, point, third, first);
  return inside || on_edge;
}

/// @brief Check whether three consecutive vertices form a convex clockwise corner.
VECGEOM_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
bool IsConvexCorner(PlanarPolygon const &polygon, size_t previous, size_t current, size_t next)
{
  const auto *x = polygon.GetVertices().x();
  const auto *y = polygon.GetVertices().y();
  const Precision cross =
      (x[current] - x[previous]) * (y[next] - y[current]) - (x[next] - x[current]) * (y[current] - y[previous]);
  return cross < 0.;
}

} // namespace PlanarPolygonUtilities
} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif
