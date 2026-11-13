#include "VecGeom/volumes/utilities/BoundingEnvelope.h"

#include "VecGeom/base/Assert.h"

#include <cmath>
#include <utility>

namespace vecgeom {

namespace {

using BoundingEnvelopeConvexPiece = clipped_extent_detail::BoundingEnvelopeConvexPiece;
using BoundingEnvelopeFace        = clipped_extent_detail::BoundingEnvelopeFace;

/// Plane inequality `normal.Dot(point) <= offset`.
struct Halfspace {
  Vector3D<Precision> normal;
  Precision offset = 0.;

  bool Contains(Vector3D<Precision> const &point, Precision tolerance) const
  {
    return normal.Dot(point) <= offset + tolerance;
  }
};

/// Solve the intersection point of three clipping planes.
bool SolvePlaneTriple(Halfspace const &a, Halfspace const &b, Halfspace const &c, Vector3D<Precision> &point)
{
  const Vector3D<Precision> bc = b.normal.Cross(c.normal);
  const Precision det          = a.normal.Dot(bc);
  if (Abs(det) <= kTolerance) return false;
  point = (a.offset * bc + b.offset * c.normal.Cross(a.normal) + c.offset * a.normal.Cross(b.normal)) / det;
  return true;
}

/// Return the arithmetic centroid used only to orient convex face planes.
Vector3D<Precision> Centroid(std::vector<Vector3D<Precision>> const &vertices)
{
  Vector3D<Precision> center(0.);
  for (auto const &vertex : vertices) {
    center += vertex;
  }
  return center / static_cast<Precision>(vertices.size());
}

/// Build one inward-facing halfspace from a convex face.
bool MakeFacePlane(std::vector<Vector3D<Precision>> const &vertices, BoundingEnvelopeFace const &face,
                   Vector3D<Precision> const &inside_point, Halfspace &plane)
{
  if (face.vertices.size() < 3) return false;
  const Vector3D<Precision> origin = vertices[face.vertices[0]];

  for (size_t i = 1; i + 1 < face.vertices.size(); ++i) {
    Vector3D<Precision> normal = (vertices[face.vertices[i]] - origin).Cross(vertices[face.vertices[i + 1]] - origin);
    const Precision normal_mag = normal.Mag();
    if (normal_mag <= kTolerance) continue;
    normal /= normal_mag;
    Precision offset = normal.Dot(origin);
    if (normal.Dot(inside_point) > offset) {
      normal = -normal;
      offset = -offset;
    }
    plane.normal = normal;
    plane.offset = offset;
    return true;
  }
  return false;
}

/// Convert all non-degenerate convex-piece faces to clipping halfspaces.
std::vector<Halfspace> MakeHalfspaces(std::vector<Vector3D<Precision>> const &vertices,
                                      std::vector<BoundingEnvelopeFace> const &faces)
{
  std::vector<Halfspace> planes;
  planes.reserve(faces.size());
  const Vector3D<Precision> inside_point = Centroid(vertices);
  for (auto const &face : faces) {
    Halfspace plane;
    if (MakeFacePlane(vertices, face, inside_point, plane)) planes.push_back(plane);
  }
  return planes;
}

/// Add active voxel-limit planes to the convex-piece clipping set.
void AddLimitHalfspaces(VoxelLimits const &limits, std::vector<Halfspace> &planes)
{
  for (int axis = 0; axis < 3; ++axis) {
    const Axis clipped_axis = AxisFromIndex(axis);
    if (!limits.IsLimited(clipped_axis)) continue;
    Vector3D<Precision> upper_normal(0.);
    upper_normal[axis] = 1.;
    planes.push_back({upper_normal, limits.Max(clipped_axis)});
    Vector3D<Precision> lower_normal(0.);
    lower_normal[axis] = -1.;
    planes.push_back({lower_normal, -limits.Min(clipped_axis)});
  }
}

/// Test a point against every halfspace with a shared tolerance.
bool InsideAll(Vector3D<Precision> const &point, std::vector<Halfspace> const &planes, Precision tolerance)
{
  for (auto const &plane : planes) {
    if (!plane.Contains(point, tolerance)) return false;
  }
  return true;
}

ClippedBox ClipConvexPiece(std::vector<Vector3D<Precision>> const &vertices,
                           std::vector<BoundingEnvelopeFace> const &faces, VoxelLimits const &limits,
                           Precision tolerance)
{
  if (limits.IsEmpty()) return {};

  std::vector<Halfspace> planes = MakeHalfspaces(vertices, faces);
  if (planes.empty()) return {};
  const size_t envelope_plane_count = planes.size();
  AddLimitHalfspaces(limits, planes);
  const Precision inclusion_tolerance = Max(tolerance, kTolerance);

  ClippedBox result;
  for (auto const &vertex : vertices) {
    if (InsideAll(vertex, planes, inclusion_tolerance)) result.Include(vertex);
  }

  // Enumerating boundary-plane triples gives the exact clipped convex envelope
  // vertices for the halfspace model used by this helper.
  for (size_t i = 0; i < planes.size(); ++i) {
    for (size_t j = i + 1; j < planes.size(); ++j) {
      for (size_t k = j + 1; k < planes.size(); ++k) {
        Vector3D<Precision> point;
        if (!SolvePlaneTriple(planes[i], planes[j], planes[k], point)) continue;
        if (InsideAll(point, planes, inclusion_tolerance)) result.Include(point);
      }
    }
  }

  // A completely unlimited query over a bounded envelope is already covered by
  // envelope vertices; keep that path robust even for duplicate or degenerate
  // input faces that suppress plane triples.
  if (!result.exists && envelope_plane_count == planes.size()) {
    for (auto const &vertex : vertices) {
      result.Include(vertex);
    }
  }

  result.Inflate(tolerance);
  IntersectBoxWithLimits(result, limits);
  return result;
}

} // namespace

BoundingEnvelope::BoundingEnvelope(Vector3D<Precision> const &local_min, Vector3D<Precision> const &local_max)
{
  BoundingEnvelopeConvexPiece piece;
  piece.local_vertices = {{local_min.x(), local_min.y(), local_min.z()}, {local_max.x(), local_min.y(), local_min.z()},
                          {local_max.x(), local_max.y(), local_min.z()}, {local_min.x(), local_max.y(), local_min.z()},
                          {local_min.x(), local_min.y(), local_max.z()}, {local_max.x(), local_min.y(), local_max.z()},
                          {local_max.x(), local_max.y(), local_max.z()}, {local_min.x(), local_max.y(), local_max.z()}};
  piece.faces = {{{0, 1, 2, 3}}, {{4, 7, 6, 5}}, {{0, 4, 5, 1}}, {{1, 5, 6, 2}}, {{2, 6, 7, 3}}, {{3, 7, 4, 0}}};
  fPieces.push_back(std::move(piece));
}

void BoundingEnvelope::AddLocalPrism(const Vector3D<Precision> *lower_base, int lower_count,
                                     const Vector3D<Precision> *upper_base, int upper_count)
{
  VECGEOM_ASSERT(lower_base != nullptr && upper_base != nullptr);
  VECGEOM_ASSERT(lower_count >= 3 && lower_count == upper_count);
  if (lower_base == nullptr || upper_base == nullptr || lower_count < 3 || lower_count != upper_count) return;

  BoundingEnvelopeConvexPiece piece;
  piece.local_vertices.reserve(2 * lower_count);
  for (int i = 0; i < lower_count; ++i)
    piece.local_vertices.push_back(lower_base[i]);
  for (int i = 0; i < upper_count; ++i)
    piece.local_vertices.push_back(upper_base[i]);

  BoundingEnvelopeFace lower_face;
  BoundingEnvelopeFace upper_face;
  lower_face.vertices.reserve(lower_count);
  upper_face.vertices.reserve(upper_count);
  for (int i = 0; i < lower_count; ++i) {
    lower_face.vertices.push_back(i);
    upper_face.vertices.push_back(lower_count + i);
  }
  piece.faces.push_back(std::move(lower_face));
  piece.faces.push_back(std::move(upper_face));

  for (int i = 0; i < lower_count; ++i) {
    const int next = (i + 1) % lower_count;
    piece.faces.push_back({{i, next, lower_count + next, lower_count + i}});
  }

  fPieces.push_back(std::move(piece));
}

#ifndef VECCORE_CUDA
void BoundingEnvelope::AddLocalConicalFrustum(Precision rmin1, Precision rmax1, Precision rmin2, Precision rmax2,
                                              Precision z1, Precision z2, Precision phi_start, Precision phi_delta,
                                              int full_phi_segments)
{
  VECGEOM_ASSERT(full_phi_segments > 0);
  VECGEOM_ASSERT(phi_delta > 0.);
  VECGEOM_ASSERT(z1 != z2);
  if (full_phi_segments <= 0 || phi_delta <= 0. || z1 == z2) return;

  phi_delta      = Min(phi_delta, kTwoPi);
  int cell_count = static_cast<int>(std::ceil(phi_delta / (kTwoPi / full_phi_segments)));
  if (cell_count < 1) cell_count = 1;

  // One convex cell uses lower inner/outer, upper inner/outer vertices and the
  // same hexahedral topology as an angular frustum wedge.
  constexpr int face_counts[]  = {4, 4, 4, 4, 4, 4};
  constexpr int face_indices[] = {
      0, 1, 2, 3, // lower z face
      4, 7, 6, 5, // upper z face
      0, 4, 7, 3, // inner or collapsed-axis face
      1, 5, 6, 2, // outer circumscribed face
      0, 4, 5, 1, // phi_i face
      3, 2, 6, 7  // phi_{i+1} face
  };
  constexpr int face_index_count = sizeof(face_indices) / sizeof(face_indices[0]);

  const Precision cell_delta  = phi_delta / cell_count;
  const Precision outer_scale = 1. / std::cos(0.5 * cell_delta);
  for (int cell = 0; cell < cell_count; ++cell) {
    const Precision phi0 = phi_start + cell * cell_delta;
    const Precision phi1 = phi0 + cell_delta;
    const Vector3D<Precision> dir0(std::cos(phi0), std::sin(phi0), 0.);
    const Vector3D<Precision> dir1(std::cos(phi1), std::sin(phi1), 0.);

    Vector3D<Precision> vertices[] = {rmin1 * dir0 + Vector3D<Precision>(0., 0., z1),
                                      outer_scale * rmax1 * dir0 + Vector3D<Precision>(0., 0., z1),
                                      outer_scale * rmax1 * dir1 + Vector3D<Precision>(0., 0., z1),
                                      rmin1 * dir1 + Vector3D<Precision>(0., 0., z1),
                                      rmin2 * dir0 + Vector3D<Precision>(0., 0., z2),
                                      outer_scale * rmax2 * dir0 + Vector3D<Precision>(0., 0., z2),
                                      outer_scale * rmax2 * dir1 + Vector3D<Precision>(0., 0., z2),
                                      rmin2 * dir1 + Vector3D<Precision>(0., 0., z2)};
    AddLocalConvexPolyhedron(vertices, 8, face_indices, face_index_count, face_counts, 6);
  }
}

void BoundingEnvelope::AddLocalSphere(Precision radius, int theta_segments, int phi_segments)
{
  VECGEOM_ASSERT(radius > 0.);
  VECGEOM_ASSERT(theta_segments >= 2);
  VECGEOM_ASSERT(phi_segments >= 3);
  if (radius <= 0. || theta_segments < 2 || phi_segments < 3) return;

  const Precision theta_step = kPi / theta_segments;
  const Precision phi_step   = kTwoPi / phi_segments;
  const Precision rtheta     = radius / std::cos(0.5 * theta_step);
  const Precision rphi       = rtheta / std::cos(0.5 * phi_step);

  std::vector<Vector3D<Precision>> phi_directions;
  phi_directions.reserve(phi_segments);
  for (int phi = 0; phi < phi_segments; ++phi) {
    const Precision angle = (phi + 0.5) * phi_step;
    phi_directions.emplace_back(std::cos(angle), std::sin(angle), 0.);
  }

  auto set_ring = [&](int theta_index, std::vector<Vector3D<Precision>> &ring) {
    const Precision theta = (theta_index + 0.5) * theta_step;
    const Precision rho   = rphi * std::sin(theta);
    const Precision z     = rtheta * std::cos(theta);
    for (int phi = 0; phi < phi_segments; ++phi) {
      ring[phi] = rho * phi_directions[phi] + Vector3D<Precision>(0., 0., z);
    }
  };

  std::vector<Vector3D<Precision>> lower_ring(phi_segments);
  std::vector<Vector3D<Precision>> upper_ring(phi_segments);
  set_ring(0, lower_ring);
  for (int theta = 1; theta < theta_segments; ++theta) {
    set_ring(theta, upper_ring);
    AddLocalPrism(lower_ring.data(), phi_segments, upper_ring.data(), phi_segments);
    lower_ring.swap(upper_ring);
  }
}

void BoundingEnvelope::AddLocalEllipticalFrustum(Precision x_radius1, Precision y_radius1, Precision x_radius2,
                                                 Precision y_radius2, Precision z1, Precision z2, int full_phi_segments)
{
  VECGEOM_ASSERT(x_radius1 > 0. && y_radius1 > 0. && x_radius2 > 0. && y_radius2 > 0.);
  VECGEOM_ASSERT(full_phi_segments >= 3);
  VECGEOM_ASSERT(z1 != z2);
  if (x_radius1 <= 0. || y_radius1 <= 0. || x_radius2 <= 0. || y_radius2 <= 0. || full_phi_segments < 3 || z1 == z2)
    return;

  constexpr int face_counts[]  = {4, 4, 4, 4, 4, 4};
  constexpr int face_indices[] = {
      0, 1, 2, 3, // lower z face
      4, 7, 6, 5, // upper z face
      0, 4, 7, 3, // collapsed-axis face
      1, 5, 6, 2, // outer circumscribed face
      0, 4, 5, 1, // phi_i face
      3, 2, 6, 7  // phi_{i+1} face
  };
  constexpr int face_index_count = sizeof(face_indices) / sizeof(face_indices[0]);

  const Precision cell_delta  = kTwoPi / full_phi_segments;
  const Precision outer_scale = 1. / std::cos(0.5 * cell_delta);
  for (int cell = 0; cell < full_phi_segments; ++cell) {
    const Precision phi0 = cell * cell_delta;
    const Precision phi1 = phi0 + cell_delta;
    const Precision cos0 = std::cos(phi0);
    const Precision sin0 = std::sin(phi0);
    const Precision cos1 = std::cos(phi1);
    const Precision sin1 = std::sin(phi1);

    Vector3D<Precision> vertices[] = {
        Vector3D<Precision>(0., 0., z1),
        Vector3D<Precision>(outer_scale * x_radius1 * cos0, outer_scale * y_radius1 * sin0, z1),
        Vector3D<Precision>(outer_scale * x_radius1 * cos1, outer_scale * y_radius1 * sin1, z1),
        Vector3D<Precision>(0., 0., z1),
        Vector3D<Precision>(0., 0., z2),
        Vector3D<Precision>(outer_scale * x_radius2 * cos0, outer_scale * y_radius2 * sin0, z2),
        Vector3D<Precision>(outer_scale * x_radius2 * cos1, outer_scale * y_radius2 * sin1, z2),
        Vector3D<Precision>(0., 0., z2)};
    AddLocalConvexPolyhedron(vertices, 8, face_indices, face_index_count, face_counts, 6);
  }
}
#endif

void BoundingEnvelope::AddLocalConvexPolyhedron(const Vector3D<Precision> *vertices, int vertex_count,
                                                const int *face_indices, int face_index_count,
                                                const int *face_vertex_counts, int face_count)
{
  VECGEOM_ASSERT(vertices != nullptr && face_indices != nullptr && face_vertex_counts != nullptr);
  VECGEOM_ASSERT(vertex_count >= 4 && face_count >= 4 && face_index_count >= 3 * face_count);
  if (vertices == nullptr || face_indices == nullptr || face_vertex_counts == nullptr || vertex_count < 4 ||
      face_count < 4 || face_index_count < 3 * face_count)
    return;

  int total_face_indices = 0;
  for (int face = 0; face < face_count; ++face) {
    const int vertex_count_on_face = face_vertex_counts[face];
    VECGEOM_ASSERT(vertex_count_on_face >= 3);
    if (vertex_count_on_face < 3) return;
    VECGEOM_ASSERT(total_face_indices + vertex_count_on_face <= face_index_count);
    if (total_face_indices + vertex_count_on_face > face_index_count) return;
    for (int i = 0; i < vertex_count_on_face; ++i) {
      const int vertex_index = face_indices[total_face_indices + i];
      VECGEOM_ASSERT(vertex_index >= 0 && vertex_index < vertex_count);
      if (vertex_index < 0 || vertex_index >= vertex_count) return;
    }
    total_face_indices += vertex_count_on_face;
  }
  VECGEOM_ASSERT(total_face_indices == face_index_count);
  if (total_face_indices != face_index_count) return;

  BoundingEnvelopeConvexPiece piece;
  piece.local_vertices.reserve(vertex_count);
  for (int i = 0; i < vertex_count; ++i)
    piece.local_vertices.push_back(vertices[i]);

  piece.faces.reserve(face_count);
  int face_index_offset = 0;
  for (int face = 0; face < face_count; ++face) {
    BoundingEnvelopeFace envelope_face;
    envelope_face.vertices.reserve(face_vertex_counts[face]);
    for (int i = 0; i < face_vertex_counts[face]; ++i)
      envelope_face.vertices.push_back(face_indices[face_index_offset + i]);
    piece.faces.push_back(std::move(envelope_face));
    face_index_offset += face_vertex_counts[face];
  }

  fPieces.push_back(std::move(piece));
}

bool BoundingEnvelope::CalculateExtent(Axis axis, VoxelLimits const &limits, Transformation3D const &transform,
                                       Precision &min, Precision &max, Precision tolerance) const
{
  const ClippedBox box = CalculateClippedBox(limits, transform, tolerance);
  if (!box.exists) return false;
  const int index = AxisIndex(axis);
  min             = box.min[index];
  max             = box.max[index];
  return true;
}

ClippedBox BoundingEnvelope::CalculateClippedBox(VoxelLimits const &limits, Transformation3D const &transform,
                                                 Precision tolerance) const
{
  ClippedBox combined;
  for (auto const &piece : fPieces) {
    std::vector<Vector3D<Precision>> mother_vertices;
    mother_vertices.reserve(piece.local_vertices.size());
    for (auto const &local_vertex : piece.local_vertices) {
      mother_vertices.push_back(transform.InverseTransform(local_vertex));
    }
    combined.Include(ClipConvexPiece(mother_vertices, piece.faces, limits, tolerance));
  }
  return combined;
}

} // namespace vecgeom
