/// @file Triangulation.h
/// @brief Generic polygon triangulation helpers.

#ifndef VECGEOM_BASE_TRIANGULATION_H
#define VECGEOM_BASE_TRIANGULATION_H

#include "VecGeom/base/Config.h"
#include "VecGeom/base/Vector.h"

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

/// @brief Triangulate a simple polygon by ear clipping.
/// @details The loop clips one ear, erases its center vertex, then restarts from
/// candidate `(0,1,2)`. Callers provide the geometric predicates so existing
/// orientation and tolerance conventions remain local to the polygon type.
/// @return True when at least one triangle was produced and every ear was found.
template <typename TriangleList, typename MakeTriangle, typename IsConvexCorner, typename ContainsPoint>
bool TriangulateEarClipping(size_t nvertices, TriangleList &triangles, MakeTriangle make_triangle,
                            IsConvexCorner is_convex_corner, ContainsPoint contains_point)
{
  if (nvertices < 3) return false;

  const size_t initial_size = triangles.size();
  Vector<size_t> indices;
  for (size_t index = 0; index < nvertices; ++index)
    indices.push_back(index);

  size_t previous = 0;
  size_t current  = 1;
  size_t next     = 2;
  while (indices.size() > 2) {
    bool found_ear = false;
    for (size_t attempt = 0; attempt < indices.size(); ++attempt) {
      const size_t first  = indices[previous];
      const size_t second = indices[current];
      const size_t third  = indices[next];
      bool good_ear       = is_convex_corner(first, second, third);

      for (size_t index : indices) {
        if (!good_ear) break;
        if (index == first || index == second || index == third) continue;
        good_ear = !contains_point(first, second, third, index);
      }

      if (good_ear) {
        triangles.push_back(make_triangle(first, second, third));
        indices.erase(indices.begin() + current);
        previous  = 0;
        current   = 1;
        next      = 2;
        found_ear = true;
        break;
      }
      previous = (previous + 1) % indices.size();
      current  = (current + 1) % indices.size();
      next     = (next + 1) % indices.size();
    }
    if (!found_ear) return false;
  }
  return triangles.size() != initial_size;
}

} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif
