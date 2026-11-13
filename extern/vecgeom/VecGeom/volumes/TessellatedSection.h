// This file is part of VecGeom and is distributed under the
// conditions in the file LICENSE.txt in the top directory.
// For the full list of authors see CONTRIBUTORS.txt and `git log`.

/// \brief Declaration of a tessellated section.
/// \file volumes/TessellatedStruct.h
/// \author Mihaela Gheata (CERN/ISS)

#ifndef VECGEOM_VOLUMES_TESSELLATEDSECTION_H_
#define VECGEOM_VOLUMES_TESSELLATEDSECTION_H_

#include <VecCore/VecCore>

#include <vector>

#include "VecGeom/base/Vector.h"
#include "VecGeom/volumes/Tile.h"

namespace vecgeom {

VECGEOM_DEVICE_DECLARE_CONV_TEMPLATE(class, TessellatedSection, typename);

inline namespace VECGEOM_IMPL_NAMESPACE {

/// Navigation helper for adjacent quadrilateral facets forming a closed
/// convex section in between 2 Z planes. A base class defines the scalar navigation
/// interfaces.
template <typename T>
class TessellatedSection {
  // Here we should be able to use vecgeom::Vector
  template <typename U>
  /// Generic vector type
  using vector_t = vecgeom::Vector<U>;
  /// Quadrilateral facets for the section
  using Facet_t = QuadrilateralFacet<T>;

private:
  size_t fNfacets = 0;    ///< Number of triangle facets on the section
  T fZ            = 0;    ///< Z position of the section
  T fDz           = 0;    ///< Half-length in Z
  T fCubicVolume  = 0;    ///< Cubic volume
  T fSurfaceArea  = 0;    ///< Surface area
  Vector3D<T> fMinExtent; ///< Minimum extent
  Vector3D<T> fMaxExtent; ///< Maximum extent
  bool fSameZ = false;    ///< All facets are at same Z
  T fUpNorm   = 0;        ///< Up normal in case of sameZ (+1 or -1)

  vector_t<Facet_t> fFacets; ///< Contiguous quadrilateral facets for host/device navigation

#if defined(VECGEOM_VC) && !defined(VECCORE_CUDA)
  using HostReal_v                     = vecgeom::VectorBackend::Real_v;
  using HostBool_v                     = vecCore::Mask<HostReal_v>;
  static constexpr size_t kHostVecSize = vecCore::VectorSize<HostReal_v>();

  struct FacetBatch {
    Vector3D<HostReal_v> fNormals;        ///< Facet normals packed in host SIMD lanes
    HostReal_v fDistances;                ///< Plane distances packed in host SIMD lanes
    Vector3D<HostReal_v> fSideVectors[4]; ///< Edge side vectors packed in host SIMD lanes
    Vector3D<HostReal_v> fVertices[4];    ///< Facet vertices packed in host SIMD lanes
    size_t fIndices[kHostVecSize] = {};   ///< Original scalar facet index for each active lane
    HostBool_v fActive;                   ///< Active lanes; tail lanes are ignored by reductions
  };

  // Scalar facets are the copied runtime representation. Host Vc builds keep a
  // prepacked cache over the same facets to avoid per-query packing; CUDA uses
  // fFacets directly and therefore keeps a compact scalar device layout.
  std::vector<FacetBatch> fFacetBatches;

  VECGEOM_FORCE_INLINE
  void RebuildHostSimdCache()
  {
    fFacetBatches.clear();
    const size_t count = fFacets.size();
    if (count == 0) return;

    const size_t nbatches = (count + kHostVecSize - 1) / kHostVecSize;
    fFacetBatches.reserve(nbatches);
    const size_t last = count - 1;

    for (size_t begin = 0; begin < count; begin += kHostVecSize) {
      FacetBatch batch;
      HostReal_v activeValue = HostReal_v(0.);

      for (size_t lane = 0; lane < kHostVecSize; ++lane) {
        const size_t index    = begin + lane;
        const bool activeLane = index < count;
        Facet_t const &facet  = fFacets[activeLane ? index : last];

        vecCore::Set(batch.fNormals.x(), lane, facet.fNormal.x());
        vecCore::Set(batch.fNormals.y(), lane, facet.fNormal.y());
        vecCore::Set(batch.fNormals.z(), lane, facet.fNormal.z());
        vecCore::Set(batch.fDistances, lane, facet.fDistance);
        vecCore::Set(activeValue, lane, activeLane ? T(1.) : T(0.));
        batch.fIndices[lane] = activeLane ? index : count;

        for (size_t ivert = 0; ivert < 4; ++ivert) {
          vecCore::Set(batch.fSideVectors[ivert].x(), lane, facet.fSideVectors[ivert].x());
          vecCore::Set(batch.fSideVectors[ivert].y(), lane, facet.fSideVectors[ivert].y());
          vecCore::Set(batch.fSideVectors[ivert].z(), lane, facet.fSideVectors[ivert].z());
          vecCore::Set(batch.fVertices[ivert].x(), lane, facet.fVertices[ivert].x());
          vecCore::Set(batch.fVertices[ivert].y(), lane, facet.fVertices[ivert].y());
          vecCore::Set(batch.fVertices[ivert].z(), lane, facet.fVertices[ivert].z());
        }
      }
      batch.fActive = activeValue > HostReal_v(0.5);
      fFacetBatches.push_back(batch);
    }
  }

  VECGEOM_FORCE_INLINE
  static HostReal_v DistPlanes(FacetBatch const &batch, Vector3D<HostReal_v> const &point)
  {
    return point.Dot(batch.fNormals) + batch.fDistances;
  }

  VECGEOM_FORCE_INLINE
  static void InsideBatch(FacetBatch const &batch, Vector3D<HostReal_v> const &point, HostBool_v &inside)
  {
    inside = batch.fActive;
    for (size_t ivert = 0; ivert < 4; ++ivert) {
      const HostReal_v saf = (point - batch.fVertices[ivert]).Dot(batch.fSideVectors[ivert]);
      inside &= saf > HostReal_v(-kTolerance);
    }
  }

  VECGEOM_FORCE_INLINE
  static bool DistanceToInConvexBatch(FacetBatch const &batch, Vector3D<HostReal_v> const &point,
                                      Vector3D<HostReal_v> const &direction, T &distance, T &limit)
  {
    const HostReal_v proj       = NonZero(direction.Dot(batch.fNormals));
    const HostBool_v movingAway = proj > HostReal_v(-kTolerance);
    const HostBool_v movingOut  = proj > HostReal_v(kTolerance);
    const HostReal_v pdist      = DistPlanes(batch, point);
    const HostBool_v outside    = pdist > HostReal_v(-kTolerance);

    if (!vecCore::MaskEmpty(batch.fActive && outside && movingAway)) return false;

    // Entry/exit interval bounds are reduced across packed side facets.
    const HostBool_v fromOutside = batch.fActive && outside && !movingAway;
    const HostBool_v fromInside  = batch.fActive && !outside && movingOut;
    const HostReal_v dist        = -pdist / proj;

    HostReal_v dmin = -InfinityLength<HostReal_v>();
    HostReal_v dmax = InfinityLength<HostReal_v>();
    vecCore__MaskedAssignFunc(dmin, fromOutside, dist);
    vecCore__MaskedAssignFunc(dmax, fromInside, dist);
    distance = Max(distance, vecCore::ReduceMax(dmin));
    limit    = Min(limit, vecCore::ReduceMin(dmax));
    return true;
  }

  VECGEOM_FORCE_INLINE
  static bool DistanceToOutConvexBatch(FacetBatch const &batch, Vector3D<HostReal_v> const &point,
                                       Vector3D<HostReal_v> const &direction, T &distance)
  {
    const HostReal_v saf   = DistPlanes(batch, point);
    const HostBool_v valid = !batch.fActive || saf < HostReal_v(kTolerance);
    if (!vecCore::MaskFull(valid)) return false;

    const HostReal_v proj      = direction.Dot(batch.fNormals);
    const HostBool_v candidate = batch.fActive && proj > HostReal_v(0.);
    HostReal_v dist            = InfinityLength<HostReal_v>();
    vecCore__MaskedAssignFunc(dist, candidate, -saf / NonZero(proj));
    distance = Min(distance, vecCore::ReduceMin(dist));
    return true;
  }

  template <bool RecordFacet>
  VECGEOM_FORCE_INLINE static bool DistanceToInOutBatch(FacetBatch const &batch, Vector3D<HostReal_v> const &point,
                                                        Vector3D<HostReal_v> const &direction, T &dtoin, T &dtoout,
                                                        int &exitFacet)
  {
    const HostReal_v projdir = NonZero(direction.Dot(batch.fNormals));
    // Ignore roundoff-sized projections: dividing a tolerated plane offset by
    // one can manufacture a large backward distance for a tangential ray.
    const HostBool_v moving  = projdir > HostReal_v(kTolerance);
    const HostReal_v pdist   = DistPlanes(batch, point);
    const HostBool_v outside = pdist > HostReal_v(kTolerance);
    if (!vecCore::MaskEmpty(batch.fActive && outside && moving)) return false;

    // Track the ray interval through all lateral halfspaces, packed by facet.
    const HostBool_v fromInside  = batch.fActive && !outside && moving;
    const HostBool_v fromOutside = batch.fActive && outside && !moving;
    const HostReal_v dist        = -pdist / projdir;
    HostReal_v dtoin_v           = -InfinityLength<HostReal_v>();
    HostReal_v dtoout_v          = InfinityLength<HostReal_v>();
    vecCore__MaskedAssignFunc(dtoin_v, fromOutside, dist);
    vecCore__MaskedAssignFunc(dtoout_v, fromInside, dist);
    dtoin = Max(dtoin, vecCore::ReduceMax(dtoin_v));
    if constexpr (RecordFacet) {
      // Keep the packed reduction; extract a lane only for a winning batch.
      const T batchDistance = dtoout_v.min();
      if (batchDistance < dtoout) {
        dtoout    = batchDistance;
        exitFacet = batch.fIndices[(dtoout_v == HostReal_v(batchDistance)).firstOne()];
      }
    } else {
      dtoout = Min(dtoout, vecCore::ReduceMin(dtoout_v));
    }
    return true;
  }

  VECGEOM_FORCE_INLINE
  static HostReal_v DistanceToLineSegmentSquaredMasked(Vector3D<HostReal_v> const &corner0,
                                                       Vector3D<HostReal_v> const &corner1,
                                                       Vector3D<HostReal_v> const &point, HostBool_v const &mask)
  {
    HostReal_v result = InfinityLength<HostReal_v>();

    Vector3D<HostReal_v> line       = corner1 - corner0;
    Vector3D<HostReal_v> projection = point - corner0;
    HostReal_v dot0                 = projection.Dot(line);
    HostBool_v condition            = dot0 <= HostReal_v(0.) && mask;
    vecCore__MaskedAssignFunc(result, condition, (point - corner0).Mag2());
    if (vecCore::MaskFull(condition || !mask)) return result;

    HostReal_v dot1 = line.Mag2();
    condition       = dot1 <= dot0 && mask;
    vecCore__MaskedAssignFunc(result, condition, (point - corner1).Mag2());
    condition = result < InfinityLength<HostReal_v>();
    if (vecCore::MaskFull(condition || !mask)) return result;

    vecCore__MaskedAssignFunc(result, !condition && mask, ((corner0 + (dot0 / dot1) * line) - point).Mag2());
    return result;
  }

  template <bool ToIn>
  VECGEOM_FORCE_INLINE static void SafetySqBatch(FacetBatch const &batch, Vector3D<HostReal_v> const &point,
                                                 T &safetysq, int &isurf)
  {
    HostReal_v safetyv = DistPlanes(batch, point);

    Vector3D<HostReal_v> intersectionv = point - safetyv * batch.fNormals;
    HostBool_v withinBound;
    InsideBatch(batch, intersectionv, withinBound);
    if (ToIn)
      withinBound &= safetyv > HostReal_v(-kTolerance);
    else
      withinBound &= safetyv < HostReal_v(kTolerance);

    safetyv *= safetyv;
    const HostBool_v outsideFacet        = batch.fActive && !withinBound;
    Vector3D<HostReal_v> safetyvOutbound = InfinityLength<HostReal_v>();
    for (size_t ivert = 0; ivert < 4; ++ivert) {
      safetyvOutbound[ivert] = DistanceToLineSegmentSquaredMasked(
          batch.fVertices[ivert], batch.fVertices[(ivert + 1) % 4], point, outsideFacet);
    }
    vecCore__MaskedAssignFunc(safetyv, outsideFacet, safetyvOutbound.Min());

    for (size_t lane = 0; lane < kHostVecSize; ++lane) {
      if (!vecCore::Get(batch.fActive, lane)) continue;
      const T saflane = vecCore::Get(safetyv, lane);
      if (saflane < safetysq) {
        safetysq = saflane;
        isurf    = batch.fIndices[lane];
      }
    }
  }
#endif

protected:
  /// Method adding a facet to the structure.

  /** The vertices are added to the list of all vertices (including duplications)
  and the extent is re-adjusted.*/
  /// @param facet Facet to be added
  VECCORE_ATT_HOST_DEVICE
  void AddFacet(Facet_t const &facet)
  {
    // Method adding a facet to the structure. The vertices are added to the
    // list of all vertices (including duplications) and the extent is re-adjusted.
    if (fSameZ) {
      Vector3D<T> const &normal = facet.GetNormal();
      VECGEOM_ASSERT(normal.Perp() < kTolerance);
      if (fUpNorm == 0.) fUpNorm = CopySign(T(1.), normal.z());
      VECGEOM_ASSERT(fUpNorm * normal.z() > 0);
    }
    fFacets.push_back(facet);
    // Adjust extent
    T xmin =
        Min(Min(facet.fVertices[0].x(), facet.fVertices[1].x()), Min(facet.fVertices[2].x(), facet.fVertices[3].x()));
    T ymin =
        Min(Min(facet.fVertices[0].y(), facet.fVertices[1].y()), Min(facet.fVertices[2].y(), facet.fVertices[3].y()));
    fMinExtent[0] = Min(fMinExtent[0], xmin);
    fMinExtent[1] = Min(fMinExtent[1], ymin);
    fMinExtent[2] = fZ - fDz;
    T xmax =
        Max(Max(facet.fVertices[0].x(), facet.fVertices[1].x()), Max(facet.fVertices[2].x(), facet.fVertices[3].x()));
    T ymax =
        Max(Max(facet.fVertices[0].y(), facet.fVertices[1].y()), Max(facet.fVertices[2].y(), facet.fVertices[3].y()));
    fMaxExtent[0]  = Max(fMaxExtent[0], xmax);
    fMaxExtent[1]  = Max(fMaxExtent[1], ymax);
    fMaxExtent[2]  = fZ + fDz;
    size_t nfacets = fFacets.size();
    VECGEOM_VALIDATE(nfacets <= fNfacets, << "Cannot add extra facets to section");
    if (nfacets == fNfacets) {
      VECGEOM_ASSERT(CalculateConvexity() == true);
#if defined(VECGEOM_VC) && !defined(VECCORE_CUDA)
      RebuildHostSimdCache();
#endif
    }
  }

  /// Calculate convexity of the section with respect to itself
  VECCORE_ATT_HOST_DEVICE
  bool CalculateConvexity()
  {
    size_t nconvex = 0;
    for (size_t i = 0; i < fNfacets; ++i) {
      Facet_t &facet = fFacets[i];
      bool convex    = true;
      for (size_t j = 0; j < fNfacets; ++j) {
        if (j == i) continue;
        for (size_t ivert = 0; ivert < 4; ++ivert) {
          convex &= facet.DistPlane(fFacets[j].fVertices[ivert]) < kTolerance;
        }
        if (!convex) continue;
      }
      facet.fConvex = convex;
      if (convex) nconvex++;
    }
    if (nconvex == fNfacets) return true;
    return false;
  }

public:
  /// Constructor:
  /// @param nfacets Number of facets in the section
  /// @param zmin Minimum Z position
  /// @param zmax Maximum Z position
  VECCORE_ATT_HOST_DEVICE
  TessellatedSection(int nfacets, T zmin, T zmax) : fNfacets(nfacets), fZ(0.5 * (zmin + zmax)), fDz(0.5 * (zmax - zmin))
  {
    VECGEOM_VALIDATE(zmax >= zmin, << "zmin is greater than zmax");
    if (fDz < kTolerance) {
      fSameZ = true;
      fDz    = 0.;
    }
    fMinExtent.Set(InfinityLength<T>());
    fMaxExtent.Set(-InfinityLength<T>());
    fFacets.reserve(nfacets);
  }

  /// Construct a section over a pre-copied facet buffer.
  /// @param facets Preallocated contiguous facet buffer
  /// @param nfacets Number of facets in the section
  /// @param zmin Minimum Z position
  /// @param zmax Maximum Z position
  VECCORE_ATT_HOST_DEVICE
  TessellatedSection(Facet_t *facets, size_t nfacets, T zmin, T zmax)
      : fNfacets(nfacets), fZ(0.5 * (zmin + zmax)), fDz(0.5 * (zmax - zmin)), fFacets(facets, static_cast<int>(nfacets))
  {
    VECGEOM_VALIDATE(zmax >= zmin, << "zmin is greater than zmax");
    if (fDz < kTolerance) {
      fSameZ = true;
      fDz    = 0.;
    }
    fMinExtent.Set(InfinityLength<T>());
    fMaxExtent.Set(-InfinityLength<T>());
    // Device copies wrap an existing facet buffer; rebuild only cheap metadata.
    for (size_t i = 0; i < nfacets; ++i) {
      Facet_t const &facet = fFacets[i];
      if (fSameZ) {
        Vector3D<T> const &normal = facet.GetNormal();
        VECGEOM_ASSERT(normal.Perp() < kTolerance);
        if (fUpNorm == 0.) fUpNorm = CopySign(T(1.), normal.z());
        VECGEOM_ASSERT(fUpNorm * normal.z() > 0);
      }
      T xmin =
          Min(Min(facet.fVertices[0].x(), facet.fVertices[1].x()), Min(facet.fVertices[2].x(), facet.fVertices[3].x()));
      T ymin =
          Min(Min(facet.fVertices[0].y(), facet.fVertices[1].y()), Min(facet.fVertices[2].y(), facet.fVertices[3].y()));
      T xmax =
          Max(Max(facet.fVertices[0].x(), facet.fVertices[1].x()), Max(facet.fVertices[2].x(), facet.fVertices[3].x()));
      T ymax =
          Max(Max(facet.fVertices[0].y(), facet.fVertices[1].y()), Max(facet.fVertices[2].y(), facet.fVertices[3].y()));
      fMinExtent[0] = Min(fMinExtent[0], xmin);
      fMinExtent[1] = Min(fMinExtent[1], ymin);
      fMaxExtent[0] = Max(fMaxExtent[0], xmax);
      fMaxExtent[1] = Max(fMaxExtent[1], ymax);
    }
    fMinExtent[2] = fZ - fDz;
    fMaxExtent[2] = fZ + fDz;
    VECGEOM_ASSERT(CalculateConvexity() == true);
#if defined(VECGEOM_VC) && !defined(VECCORE_CUDA)
    RebuildHostSimdCache();
#endif
  }

  /// Method for adding a new quadrilateral facet
  /** @param vt0      First vertex
      @param vt1      Second vertex
      @param vt2      Third vertex
      @param vt3      Fourth vertex
      @param absolute If true then vt0, vt1, vt2 and vt3 are the vertices to be added in
        anti-clockwise order looking from the outsider. If false the vertices are relative
        to the first: vt0, vt0+vt1, vt0+vt2, vt0+vt3 in anti-clockwise order when looking from the
        outsider.
  */
  VECCORE_ATT_HOST_DEVICE
  bool AddQuadrilateralFacet(Vector3D<T> const &vt0, Vector3D<T> const &vt1, Vector3D<T> const &vt2,
                             Vector3D<T> const &vt3, bool absolute = true)
  {
    // Quadrilateral facet, normal pointing outside
    Facet_t facet;
    if (absolute) {
      if (!facet.SetVertices(vt0, vt1, vt2, vt3)) {
        return false;
      }
      AddFacet(facet);
    } else {
      if (!facet.SetVertices(vt0, vt0 + vt1, vt0 + vt1 + vt2, vt0 + vt1 + vt2 + vt3)) {
        return false;
      }
      AddFacet(facet);
    }
    return true;
  }

  /// Fast check using the extent if the point is outside
  /// @param point Point position
  /// @return True if point is outside
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  bool IsOutside(Vector3D<T> const &point)
  {
    return ((point - fMinExtent).Min() < -kTolerance || (point - fMaxExtent).Max() > kTolerance);
  }
  /// Check if a point is inside the section
  /// @param point Point position
  /// @return Enumeration inside/outside/surface
  template <typename Real_v>
  VECCORE_ATT_HOST_DEVICE VECGEOM_FORCE_INLINE Inside_t Inside(Vector3D<Real_v> const &point) const
  {
#if defined(VECGEOM_VC) && !defined(VECCORE_CUDA)
    if constexpr (vecCore::VectorSize<Real_v>() == 1) {
      if (fFacetBatches.size() != 0) {
        bool inside = true;
        Vector3D<HostReal_v> pointv(point);
        for (size_t ibatch = 0; ibatch < fFacetBatches.size(); ++ibatch) {
          FacetBatch const &batch     = fFacetBatches[ibatch];
          const HostReal_v distPlanes = DistPlanes(batch, pointv);
          if (!vecCore::MaskEmpty(batch.fActive && distPlanes > HostReal_v(kTolerance))) return kOutside;
          inside &= vecCore::MaskFull(!batch.fActive || distPlanes < HostReal_v(-kTolerance));
        }
        return inside ? kInside : kSurface;
      }
    }
#endif
    // Assuming the fast check on extent was already done using the scalar point
    bool inside = true;
    for (size_t i = 0; i < fFacets.size(); ++i) {
      Real_v distPlanes = fFacets[i].DistPlane(point);
      if (distPlanes > Real_v(kTolerance)) return kOutside;
      inside &= distPlanes < Real_v(-kTolerance);
    }
    if (inside) return kInside;
    return kSurface;
  }

  /// Getter for the number of facets of the section
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  size_t GetNfacets() const { return fFacets.size(); }

  /// Getter for the facet at a given index
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Facet_t const &GetFacet(size_t i) const { return fFacets[i]; }

  /// @brief Record whether a facet's outward plane supports the complete owner solid.
  /// @details Extruded construction promotes the section-local convexity bit
  /// after all sections are available; exit-normal queries then read it without
  /// a second support array or runtime facet-index mapping.
  /// @param i Section facet index.
  /// @param supports Whole-owner supporting-plane classification.
  VECCORE_ATT_HOST_DEVICE
  void SetFacetSupportsSolid(size_t i, bool supports) { fFacets[i].fConvex = supports; }

  /// Getter for the contiguous facet buffer.
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Facet_t const *GetFacetData() const { return fFacets.size() == 0 ? nullptr : &fFacets[0]; }

  /// Check if point is inside the section. Note that the Z range is not checked.
  /// @param point Point position
  // @return True if point is inside or on the surface
  template <typename Real_v>
  VECCORE_ATT_HOST_DEVICE VECGEOM_FORCE_INLINE bool Contains(Vector3D<Real_v> const &point) const
  {
#if defined(VECGEOM_VC) && !defined(VECCORE_CUDA)
    if constexpr (vecCore::VectorSize<Real_v>() == 1) {
      if (fFacetBatches.size() != 0) {
        Vector3D<HostReal_v> pointv(point);
        for (size_t ibatch = 0; ibatch < fFacetBatches.size(); ++ibatch) {
          FacetBatch const &batch = fFacetBatches[ibatch];
          if (!vecCore::MaskEmpty(batch.fActive && DistPlanes(batch, pointv) > HostReal_v(0.))) return false;
        }
        return true;
      }
    }
#endif
    // Check if point is in the bounding box
    // if ((point - fMinExtent).Min() < 0. || (point - fMaxExtent).Max() > 0.) return kOutside;

    for (size_t i = 0; i < fFacets.size(); ++i) {
      if (fFacets[i].DistPlane(point) > Real_v(0)) return false;
    }
    return true;
  }

  /// Computes distance to in for the cluster, but may ignore Z checks
  /// @tparam knownFacets Positive when the caller knows the exact facet count,
  /// allowing the compiler to specialize the scalar loop bound.
  /// @param [in]  point Point position
  /// @param [in]  direction Direction for the distance computation
  /// @param [in]  invdirz Inverse of direction on Z
  /// @param [in]  stepmax Search limit for the distance
  /// @return Computed distance
  template <bool skipZ = true, int knownFacets = 0>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE T DistanceToIn(Vector3D<T> const &point, Vector3D<T> const &direction,
                                                              T invdirz, T stepmax) const
  {
    // Compute distance to segment from outside point.
    if (fSameZ) {
      // All facets are on the plane at fZ
      // Distance to plane
      T pz = point.z() - fZ;
      // If wrong direction or opposite side, no hit
      if (fUpNorm * direction.z() > 0 || pz * fUpNorm < -kTolerance) return InfinityLength<T>();
      T distance = -pz * invdirz;
      // Still need to check that the propagated point is in the section
      Vector3D<T> propagated(point + distance * direction);
      for (size_t i = 0; i < fFacets.size(); ++i) {
        if (fFacets[i].Contains(propagated)) return distance;
      }
      return InfinityLength<T>();
    }

    T pz = point.z() - fZ;
    if (!skipZ) {
      if ((Abs(pz) - fDz) > -kTolerance && pz * direction.z() >= 0) return InfinityLength<T>();
    }
    const T ddz   = CopySign(fDz, invdirz);
    const T distz = -(pz + ddz) * invdirz;
    T distance    = distz;
    T limit       = Min(-(pz - ddz) * invdirz, stepmax);

#if defined(VECGEOM_VC) && !defined(VECCORE_CUDA)
    if (fFacetBatches.size() != 0) {
      Vector3D<HostReal_v> pointv(point);
      Vector3D<HostReal_v> directionv(direction);
      for (size_t ibatch = 0; ibatch < fFacetBatches.size(); ++ibatch) {
        if (!DistanceToInConvexBatch(fFacetBatches[ibatch], pointv, directionv, distance, limit))
          return InfinityLength<T>();
        if (distance >= limit - kTolerance) return InfinityLength<T>();
      }
      if (skipZ) {
        if (distance > distz) return distance;
        return InfinityLength<T>();
      }
      return distance;
    }
#endif

    // Convex section entry is the largest lower ray interval bound; exit is
    // the smallest upper bound. If any side plane inverts the interval, the ray
    // misses the section.
    const size_t nfacets = knownFacets > 0 ? knownFacets : fFacets.size();
    for (size_t i = 0; i < nfacets; ++i) {
      const T proj          = NonZero(direction.Dot(fFacets[i].fNormal));
      const bool movingAway = proj > -kTolerance;
      // A tangent can reject an outside entry but cannot bound an inside exit.
      const bool movingOut = proj > kTolerance;
      const T pdist        = fFacets[i].DistPlane(point);
      const bool outside   = pdist > -kTolerance;
      if (outside && movingAway) return InfinityLength<T>();

      const T dist = -pdist / proj;
      if (outside) {
        if (dist > distance) distance = dist;
      } else if (movingOut) {
        if (dist < limit) limit = dist;
      }
    }
    // Entry can only increase and the exit limit can only decrease, so an
    // inverted interval cannot become valid again in a later iteration.
    if (distance >= limit - kTolerance) return InfinityLength<T>();
    if (skipZ) {
      if (distance > distz) return distance;
      return InfinityLength<T>();
    }
    return distance;
  }

  /// Computes distance to out for the cluster, but may ignore Z checks
  /// @param [in]  point Point position
  /// @param [in]  direction Direction for the distance computation
  /// @return Computed distance
  template <bool skipZ = true>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE T DistanceToOut(Vector3D<T> const &point,
                                                               Vector3D<T> const &direction) const
  {
    // Compute distance to segment from point inside, returning also the crossed
    // facet.
    T pz         = point.z() - fZ;
    const T safz = Abs(pz) - fDz;
    if (safz > kTolerance) return -kTolerance;
    const T vz = direction.z();
    T distance = vz == T(0.) ? InfinityLength<T>() : (CopySign(fDz, vz) - pz) / vz;
#if defined(VECGEOM_VC) && !defined(VECCORE_CUDA)
    if (fFacetBatches.size() != 0) {
      Vector3D<HostReal_v> pointv(point);
      Vector3D<HostReal_v> directionv(direction);
      for (size_t ibatch = 0; ibatch < fFacetBatches.size(); ++ibatch) {
        if (!DistanceToOutConvexBatch(fFacetBatches[ibatch], pointv, directionv, distance)) return -kTolerance;
      }
      return Max(distance, T(0.));
    }
#endif
    T dist;
    // For an inside point, every side plane must be on the accepted halfspace;
    // the exit is the nearest plane with an outward direction projection.
    for (size_t i = 0; i < fFacets.size(); ++i) {
      const T saf = fFacets[i].DistPlane(point);
      if (saf >= kTolerance) return -kTolerance;
      const T proj = direction.Dot(fFacets[i].fNormal);
      if (proj > T(0.)) {
        dist = -saf / NonZero(proj);
        if (dist < distance) distance = dist;
      }
    }
    // A tolerated outside-plane offset can have a far backward root at a
    // shallow angle. It is a surface exit, not a wrong-side distance.
    return Max(distance, T(0.));
  }

  /// Compute distance to segment from point inside, returning also the crossed facet.
  /// @param [in]  point Point position
  /// @param [in]  direction Direction for the distance computation
  /// @param [in]  invdirz Inverse of direction on Z
  /// @param[out] hitFacet Optional selected section facet; required when RecordFacet is true.
  /// @return Computed distance
  template <bool RecordFacet = false>
  VECCORE_ATT_HOST_DEVICE VECGEOM_FORCE_INLINE T DistanceToOutRange(Vector3D<T> const &point,
                                                                    Vector3D<T> const &direction, T invdirz,
                                                                    int *hitFacet = nullptr) const
  {
    int exitFacet = -1;
    if constexpr (RecordFacet) *hitFacet = -1;
    // Compute distance to segment from point inside, returning also the crossed
    // facet.
    if (fSameZ) {
      // All facets are on the plane at z = fZ
      // Distance to plane
      T pz = point.z() - fZ;
      // If wrong direction or opposite side, no hit
      if (fUpNorm * direction.z() < 0 || pz * fUpNorm > kTolerance) return InfinityLength<T>();
      T distance = -pz * invdirz;
      // Still need to check that the propagated point is in the section
      Vector3D<T> propagated(point + distance * direction);
      for (size_t i = 0; i < fFacets.size(); ++i) {
        if (fFacets[i].Contains(propagated)) {
          if constexpr (RecordFacet) *hitFacet = static_cast<int>(i);
          return distance;
        }
      }
      return InfinityLength<T>();
    }
    T pz     = point.z() - fZ;
    T dmax   = (CopySign(fDz, invdirz) - pz) * invdirz;
    T dmin   = (-CopySign(fDz, invdirz) - pz) * invdirz;
    T dtoin  = dmin;                // will be reduced Max for all clusters
    T dtoout = InfinityLength<T>(); // will be reduced Min for all clusters
#if defined(VECGEOM_VC) && !defined(VECCORE_CUDA)
    if (fFacetBatches.size() != 0) {
      Vector3D<HostReal_v> pointv(point);
      Vector3D<HostReal_v> directionv(direction);
      for (size_t ibatch = 0; ibatch < fFacetBatches.size(); ++ibatch) {
        if (!DistanceToInOutBatch<RecordFacet>(fFacetBatches[ibatch], pointv, directionv, dtoin, dtoout, exitFacet))
          return InfinityLength<T>();
      }
      if (dtoout > dtoin - kTolerance && dtoout < dmax) {
        if constexpr (RecordFacet) *hitFacet = exitFacet;
        return dtoout;
      }
      return InfinityLength<T>();
    }
#endif
    // Track the ray interval through all lateral halfspaces of the convex
    // section. Entering bounds raise dtoin, exiting bounds lower dtoout.
    for (size_t i = 0; i < fFacets.size(); ++i) {
      const T projdir = NonZero(direction.Dot(fFacets[i].fNormal));
      // Ignore roundoff-sized projections for the same reason as the packed
      // facet path above.
      const bool moving  = projdir > kTolerance;
      const T projdist   = fFacets[i].DistPlane(point);
      const bool outside = projdist > kTolerance;
      if (outside && moving) return InfinityLength<T>();
      const T dist = -projdist / projdir;
      if (outside) {
        if (dist > dtoin) dtoin = dist;
      } else if (moving) {
        if (dist < dtoout) {
          dtoout    = dist;
          exitFacet = static_cast<int>(i);
        }
      }
    }

    if (dtoout > dtoin - kTolerance && dtoout < dmax) {
      if constexpr (RecordFacet) *hitFacet = exitFacet;
      return dtoout;
    }
    return InfinityLength<T>();
  }

  /// @brief Find the nearest outward lateral-plane crossing of a convex section.
  /// @param point Start point inside or on the section
  /// @param direction Unit ray direction
  /// @param[out] isurf Selected facet index, or -1 for a miss or rejected start
  /// @return Lateral exit candidate, infinity if none, or -kTolerance for a
  /// start outside a lateral half-space. Tolerated backward crossings return zero.
  /// @details The caller must combine this candidate with the Z-cap distance;
  /// no finite-facet or Z-range test is performed here.
  template <typename Real_v>
  VECCORE_ATT_HOST_DEVICE VECGEOM_FORCE_INLINE T DistanceToOut(Vector3D<Real_v> const &point,
                                                               Vector3D<Real_v> const &direction, int &isurf) const
  {
#if defined(VECGEOM_VC) && !defined(VECCORE_CUDA)
    if (fFacetBatches.size() != 0) {
      isurf      = -1;
      T distance = InfinityLength<T>();
      Vector3D<HostReal_v> pointv(point);
      Vector3D<HostReal_v> directionv(direction);

      for (size_t ibatch = 0; ibatch < fFacetBatches.size(); ++ibatch) {
        FacetBatch const &batch = fFacetBatches[ibatch];
        const HostReal_v proj   = directionv.Dot(batch.fNormals);
        const HostReal_v saf    = DistPlanes(batch, pointv);
        if (!vecCore::MaskFull(!batch.fActive || saf < HostReal_v(kTolerance))) {
          isurf = -1;
          return -kTolerance;
        }
        const HostBool_v valid = batch.fActive && proj > HostReal_v(0.);
        HostReal_v dist        = InfinityLength<HostReal_v>();
        vecCore__MaskedAssignFunc(dist, valid, -saf / NonZero(proj));

        // Use the common SSE/AVX reduction and resolve ties to the first lane.
        const T batchDistance = dist.min();
        if (batchDistance < distance) {
          distance = batchDistance;
          isurf    = batch.fIndices[(dist == HostReal_v(batchDistance)).firstOne()];
        }
      }
      return Max(distance, T(0.));
    }
#endif
    // Compute distance to segment from point inside, returning also the crossed
    // facet.
    isurf      = -1;
    T distance = InfinityLength<T>();
    for (size_t i = 0; i < fFacets.size(); ++i) {
      const Real_v proj = direction.Dot(fFacets[i].fNormal);
      const Real_v saf  = fFacets[i].DistPlane(point);
      if (saf >= Real_v(kTolerance)) {
        isurf = -1;
        return -kTolerance;
      }
      if (proj <= Real_v(0.)) continue;
      const Real_v dist = -saf / NonZero(proj);
      // The nearest outward plane crossing lies on the finite facet when it
      // precedes the cap exit; the caller performs that final minimization.
      if (dist < distance) {
        distance = dist;
        isurf    = i;
      }
    }
    return Max(distance, T(0.));
  }

  /// Compute safety from outside point
  /// @param [in] point Point position
  /// @return Safety value
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  T SafetyToIn(Vector3D<T> const &point) const
  {
    // Compute approximate safety for the convex case
    T safety = Max(fZ - fDz - point.z(), point.z() - fZ - fDz);
#if defined(VECGEOM_VC) && !defined(VECCORE_CUDA)
    if (fFacetBatches.size() != 0) {
      Vector3D<HostReal_v> pointv(point);
      for (size_t ibatch = 0; ibatch < fFacetBatches.size(); ++ibatch) {
        HostReal_v saf = -InfinityLength<HostReal_v>();
        vecCore__MaskedAssignFunc(saf, fFacetBatches[ibatch].fActive, DistPlanes(fFacetBatches[ibatch], pointv));
        safety = Max(safety, vecCore::ReduceMax(saf));
      }
      return safety;
    }
#endif
    for (size_t i = 0; i < fFacets.size(); ++i) {
      const T saf = fFacets[i].DistPlane(point);
      if (saf > safety) safety = saf;
    }
    return safety;
  }

  /// Compute safety squared from outside point and closest surface
  /// @param [in] point Point position
  /// @param [out] isurf Closest surface
  /// @return Safety squared value
  template <typename Real_v>
  VECGEOM_FORCE_INLINE T SafetyToInSq(Vector3D<Real_v> const &point, int &isurf) const
  {
#if defined(VECGEOM_VC) && !defined(VECCORE_CUDA)
    if constexpr (vecCore::VectorSize<Real_v>() == 1) {
      if (fFacetBatches.size() != 0) {
        isurf      = -1;
        T safetysq = InfinityLength<T>();
        Vector3D<HostReal_v> pointv(point);
        for (size_t ibatch = 0; ibatch < fFacetBatches.size(); ++ibatch) {
          SafetySqBatch<true>(fFacetBatches[ibatch], pointv, safetysq, isurf);
        }
        return safetysq;
      }
    }
#endif
    // Compute safety squared to segment from point outside, returning also the crossed
    // facet.
    isurf      = -1;
    T safetysq = InfinityLength<T>();
    for (size_t i = 0; i < fFacets.size(); ++i) {
      const T safsqcrt = fFacets[i].template SafetySq<true>(point);
      if (safsqcrt < safetysq) {
        safetysq = safsqcrt;
        isurf    = i;
      }
    }
    return safetysq;
  }

  /// Compute safety from point inside
  /// @param [in] point Point position
  /// @return Safety value
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  T SafetyToOut(Vector3D<T> const &point) const
  {
    // Compute approximate safety for the convex case
    T safety = Max(fZ - fDz - point.z(), point.z() - fZ - fDz);
#if defined(VECGEOM_VC) && !defined(VECCORE_CUDA)
    if (fFacetBatches.size() != 0) {
      Vector3D<HostReal_v> pointv(point);
      for (size_t ibatch = 0; ibatch < fFacetBatches.size(); ++ibatch) {
        HostReal_v saf = -InfinityLength<HostReal_v>();
        vecCore__MaskedAssignFunc(saf, fFacetBatches[ibatch].fActive, DistPlanes(fFacetBatches[ibatch], pointv));
        safety = Max(safety, vecCore::ReduceMax(saf));
      }
      return -safety;
    }
#endif
    for (size_t i = 0; i < fFacets.size(); ++i) {
      const T saf = fFacets[i].DistPlane(point);
      if (saf > safety) safety = saf;
    }
    return -safety;
  }

  /// Compute side-only safety from an inside point.
  /// @param [in] point Point position
  /// @return Safety to the nearest lateral facet, ignoring section z caps
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  T LateralSafetyToOut(Vector3D<T> const &point) const
  {
    T safety = -InfinityLength<T>();
#if defined(VECGEOM_VC) && !defined(VECCORE_CUDA)
    if (fFacetBatches.size() != 0) {
      Vector3D<HostReal_v> pointv(point);
      for (size_t ibatch = 0; ibatch < fFacetBatches.size(); ++ibatch) {
        HostReal_v saf = -InfinityLength<HostReal_v>();
        vecCore__MaskedAssignFunc(saf, fFacetBatches[ibatch].fActive, DistPlanes(fFacetBatches[ibatch], pointv));
        safety = Max(safety, vecCore::ReduceMax(saf));
      }
      return -safety;
    }
#endif
    for (size_t i = 0; i < fFacets.size(); ++i) {
      const T saf = fFacets[i].DistPlane(point);
      if (saf > safety) safety = saf;
    }
    return -safety;
  }

  /// Compute safety squared from point inside and closest surface
  /// @param [in] point Point position
  /// @param [out] isurf Closest surface
  /// @return Safety squared value
  template <typename Real_v>
  VECGEOM_FORCE_INLINE T SafetyToOutSq(Vector3D<Real_v> const &point, int &isurf) const
  {
#if defined(VECGEOM_VC) && !defined(VECCORE_CUDA)
    if constexpr (vecCore::VectorSize<Real_v>() == 1) {
      if (fFacetBatches.size() != 0) {
        isurf      = -1;
        T safetysq = InfinityLength<T>();
        Vector3D<HostReal_v> pointv(point);
        for (size_t ibatch = 0; ibatch < fFacetBatches.size(); ++ibatch) {
          SafetySqBatch<false>(fFacetBatches[ibatch], pointv, safetysq, isurf);
        }
        return safetysq;
      }
    }
#endif
    // Compute safety squared to segment from point inside, returning also the crossed
    // facet.
    isurf      = -1;
    T safetysq = InfinityLength<T>();
    for (size_t i = 0; i < fFacets.size(); ++i) {
      const T safsqcrt = fFacets[i].template SafetySq<false>(point);
      if (safsqcrt < safetysq) {
        safetysq = safsqcrt;
        isurf    = i;
      }
    }
    return safetysq;
  }

  /// Compute normal to segment surface in given point near surface.
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  void Normal(Vector3D<T> const & /*point*/, Vector3D<T> & /*normal*/, bool & /*valid*/) const {}
};

std::ostream &operator<<(std::ostream &os, TessellatedSection<double> const &ts);

} // namespace VECGEOM_IMPL_NAMESPACE
} // end namespace vecgeom

#endif
