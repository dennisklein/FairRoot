// This file is part of VecGeom and is distributed under the
// conditions in the file LICENSE.txt in the top directory.
// For the full list of authors see CONTRIBUTORS.txt and `git log`.

/// \brief Host-side clipped extent helper for AABB and convex polyhedral envelopes.
/// \file volumes/utilities/BoundingEnvelope.h

#ifndef VECGEOM_VOLUMES_UTILITIES_BOUNDINGENVELOPE_H_
#define VECGEOM_VOLUMES_UTILITIES_BOUNDINGENVELOPE_H_

#include "VecGeom/volumes/utilities/ClippedExtent.h"

#include <vector>

namespace vecgeom {

namespace clipped_extent_detail {

/// @brief One convex envelope face described by local vertex indices.
struct BoundingEnvelopeFace {
  /// Indices into the owning convex piece vertex array.
  std::vector<int> vertices;
};

/// @brief One convex local envelope piece consumed by `BoundingEnvelope`.
struct BoundingEnvelopeConvexPiece {
  /// Local-frame vertices of the convex piece.
  std::vector<Vector3D<Precision>> local_vertices;
  /// Faces defining the convex piece boundary.
  std::vector<BoundingEnvelopeFace> faces;
};

} // namespace clipped_extent_detail

/// @brief Host-side clipped extent helper for conservative convex envelopes.
///
/// @details
/// A `BoundingEnvelope` stores one or more local-frame convex pieces. Queries
/// transform those pieces to the mother frame, clip them against `VoxelLimits`,
/// and return conservative axis-aligned bounds. Multiple pieces are combined as
/// a union, so non-convex solids can be represented by conservative convex
/// decompositions.
class BoundingEnvelope {
public:
  BoundingEnvelope() = default;

  /// @brief Construct an envelope from one local AABB.
  /// @param local_min Local lower AABB corner.
  /// @param local_max Local upper AABB corner.
  BoundingEnvelope(Vector3D<Precision> const &local_min, Vector3D<Precision> const &local_max);

  /// @brief Construct an envelope from one local AABB.
  /// @param local_min Local lower AABB corner.
  /// @param local_max Local upper AABB corner.
  /// @return Bounding envelope containing one AABB piece.
  static BoundingEnvelope FromLocalAABB(Vector3D<Precision> const &local_min, Vector3D<Precision> const &local_max)
  {
    return BoundingEnvelope(local_min, local_max);
  }

  /// @brief Add one prism with matching lower and upper polygon bases.
  /// @param lower_base Local vertices of the lower base polygon.
  /// @param lower_count Number of vertices in `lower_base`.
  /// @param upper_base Local vertices of the upper base polygon.
  /// @param upper_count Number of vertices in `upper_base`; must match `lower_count`.
  void AddLocalPrism(const Vector3D<Precision> *lower_base, int lower_count, const Vector3D<Precision> *upper_base,
                     int upper_count);

#ifndef VECCORE_CUDA
  /// @brief Add a conservative angular decomposition of one local conical frustum.
  ///
  /// @details
  /// The frustum is split into convex phi cells. Outer circular arcs are
  /// circumscribed by scaling the outer radii on each cell boundary by
  /// `1 / cos(cell_delta_phi / 2)`, matching the envelope idea used by
  /// Geant4 polycone extents. Inner radii stay unscaled, which may over-cover
  /// the hole but does not under-cover the solid. `full_phi_segments` controls
  /// the number of cells used for a full `2*pi` revolution.
  ///
  /// @param rmin1 Inner radius at `z1`.
  /// @param rmax1 Outer radius at `z1`.
  /// @param rmin2 Inner radius at `z2`.
  /// @param rmax2 Outer radius at `z2`.
  /// @param z1 Lower or first z plane of the frustum.
  /// @param z2 Upper or second z plane of the frustum.
  /// @param phi_start Start angle in local coordinates.
  /// @param phi_delta Phi opening angle in radians.
  /// @param full_phi_segments Number of angular cells used for a full circle.
  void AddLocalConicalFrustum(Precision rmin1, Precision rmax1, Precision rmin2, Precision rmax2, Precision z1,
                              Precision z2, Precision phi_start, Precision phi_delta, int full_phi_segments = 24);

  /// @brief Add a conservative circumscribed envelope for a full local sphere.
  ///
  /// @details
  /// The sphere is decomposed into theta bands with regular phi polygons.
  /// Ring vertices are placed at angular cell centers and scaled by
  /// `1 / cos(half_cell_angle)`, so the meridian and phi polygons both
  /// circumscribe the spherical surface. The default 8 theta and 16 phi
  /// segments match the `G4Orb` clipped-extent envelope.
  ///
  /// @param radius Sphere radius.
  /// @param theta_segments Number of bands in polar angle.
  /// @param phi_segments Number of angular cells in a full circle.
  void AddLocalSphere(Precision radius, int theta_segments = 8, int phi_segments = 16);

  /// @brief Add a conservative full-phi elliptical frustum envelope.
  ///
  /// @details
  /// Each z end is covered by an affine-scaled tangent polygon of the local
  /// ellipse, using the same `1 / cos(cell_delta_phi / 2)` circumscription as
  /// circular frustum envelopes. The helper covers solid elliptical frusta with
  /// no inner hole; equal lower and upper radii produce an elliptical cylinder.
  ///
  /// @param x_radius1 Semi-axis in x at `z1`.
  /// @param y_radius1 Semi-axis in y at `z1`.
  /// @param x_radius2 Semi-axis in x at `z2`.
  /// @param y_radius2 Semi-axis in y at `z2`.
  /// @param z1 Lower or first z plane of the frustum.
  /// @param z2 Upper or second z plane of the frustum.
  /// @param full_phi_segments Number of angular cells used for a full ellipse.
  void AddLocalEllipticalFrustum(Precision x_radius1, Precision y_radius1, Precision x_radius2, Precision y_radius2,
                                 Precision z1, Precision z2, int full_phi_segments = 24);
#endif

  /// @brief Add one finite convex polyhedral envelope piece from local vertices and faces.
  ///
  /// @details
  /// Faces are described by contiguous entries in `face_indices`; each
  /// `face_vertex_counts[i]` gives the number of vertex indices for face `i`.
  /// Face orientation is not significant: plane normals are oriented from the
  /// piece centroid when the envelope is clipped. This helper is intended for
  /// exact or conservative convex pieces such as tetrahedra, pyramids, and
  /// future convex decompositions; non-convex solids must add multiple convex
  /// pieces or use another conservative envelope.
  ///
  /// @param vertices Local-frame vertices.
  /// @param vertex_count Number of vertices.
  /// @param face_indices Contiguous face vertex indices.
  /// @param face_index_count Number of entries in `face_indices`.
  /// @param face_vertex_counts Number of vertices for each face.
  /// @param face_count Number of faces.
  void AddLocalConvexPolyhedron(const Vector3D<Precision> *vertices, int vertex_count, const int *face_indices,
                                int face_index_count, const int *face_vertex_counts, int face_count);

  /// @brief Calculate one clipped mother-frame extent interval.
  /// @param axis Axis to project.
  /// @param limits Mother-frame clipping limits.
  /// @param transform Local-to-mother transform.
  /// @param[out] min Lower extent along `axis` when an overlap exists.
  /// @param[out] max Upper extent along `axis` when an overlap exists.
  /// @param tolerance Final symmetric inflation applied before limiting.
  /// @return True when the envelope intersects `limits`.
  bool CalculateExtent(Axis axis, VoxelLimits const &limits, Transformation3D const &transform, Precision &min,
                       Precision &max, Precision tolerance = kAABBPadding) const;

  /// @brief Calculate a clipped mother-frame AABB for all stored pieces.
  /// @param limits Mother-frame clipping limits.
  /// @param transform Local-to-mother transform.
  /// @param tolerance Final symmetric inflation applied before limiting.
  /// @return Conservative clipped box; `exists` is false when no piece overlaps.
  ClippedBox CalculateClippedBox(VoxelLimits const &limits, Transformation3D const &transform,
                                 Precision tolerance = kAABBPadding) const;

private:
  std::vector<clipped_extent_detail::BoundingEnvelopeConvexPiece> fPieces;
};

} // namespace vecgeom

#endif // VECGEOM_VOLUMES_UTILITIES_BOUNDINGENVELOPE_H_
