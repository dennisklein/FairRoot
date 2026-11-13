// This file is part of VecGeom and is distributed under the
// conditions in the file LICENSE.txt in the top directory.
// For the full list of authors see CONTRIBUTORS.txt and `git log`.

/// \brief Declaration of the data structure for the tessellated shape.
/// \file volumes/TessellatedStruct.h
/// \author First version created by Mihaela Gheata (CERN/ISS)

#ifndef VECGEOM_VOLUMES_TESSELLATEDSTRUCT_H_
#define VECGEOM_VOLUMES_TESSELLATEDSTRUCT_H_

#include "VecGeom/base/RNG.h"
#include "VecGeom/base/Config.h"
#include "VecGeom/volumes/Tile.h"
#include <vector>

#include "VecGeom/base/Stopwatch.h"
#include "VecGeom/base/IndexedBVH.h"
#include "VecGeom/management/ABBoxManager.h"

namespace vecgeom {

// VECGEOM_DEVICE_DECLARE_CONV_TEMPLATE(class, TessellatedStruct, typename);
VECGEOM_DEVICE_DECLARE_CONV_TEMPLATE_1v_1t(class, TessellatedStruct, size_t, typename);

inline namespace VECGEOM_IMPL_NAMESPACE {

/// Compact explicit-index hierarchy shared by tessellated solids and the faceted extruded helper.
using FacetBVH_t = IndexedBVH<float>;

/** Templated class holding the data structures for the tessellated solid.

  The class is templated on the number of edges for the composing facets and on the floating precision type
  used to represent the data and perform all calculations. It provides API for:
  - Adding triangular or quadrilateral facets
  - Initializing internally all data structures after adding all the facets that compose a closed tessellated
  surface. This creates a temporary GridHelper used for vertex de-duplication and then constructs the runtime BVH.
*/
template <size_t NVERT, typename T = double>
class TessellatedStruct {

  using Facet_t = Tile<NVERT, T>;

  // Here we should be able to use vecgeom::Vector
  template <typename U>
  using vector_t = vecgeom::Vector<U>;

  /** Structure representing a cell of a uniform grid embedding the tessellated solid bounding volume.

  The cell just stores an array of vertex indices. The cell coordinates are known by the
  GridHelper data structure owning it.
  */
  struct GridCell {
    vector_t<int> fArray; ///< Array of vertex indices

    /// Default constructor for a grid cell.
    VECCORE_ATT_HOST_DEVICE
    GridCell() { /* fArray.reserve(4); */ }
  };

  /** Helper structure representing a grid of equal cells dividing the bounding box of a tessellated solid.

  The grid helper is defined by the number of cells in X/Y/Z. It provides methods to retrieve the cell
  containing a space point, or associated to a triplet of indices on (x, y, z). The helper is used by the
  TessellatedStruct while removing duplicate vertices during construction.
  */
  //__________________________________________________________________________
  struct GridHelper {
    int fNgrid       = 0;           ///< Grid size
    int fNcells      = 0;           ///< Number of cells in the grid
    GridCell **fGrid = nullptr;     ///< Grid for vertex de-duplication
    Vector3D<T> fMinExtent;         ///< Minimum extent
    Vector3D<T> fMaxExtent;         ///< Maximum extent
    Vector3D<T> fInvExtSize;        ///< Inverse extent size
    vector_t<Vector3D<T>> fAllVert; ///< Full list of vertices

    /// Default constructor for the grid helper structure.
    VECCORE_ATT_HOST_DEVICE
    GridHelper() {}

    /// Destructor of the grid helper, deleting the cells and their content.
    VECCORE_ATT_HOST_DEVICE
    ~GridHelper()
    {
      if (fGrid) {
        for (int i = 0; i < fNcells; ++i)
          delete fGrid[i];
        delete[] fGrid;
      }
    }

    /// Create all the cells corresponding to a given grid division number.
    /** @param ngrid Number of cells on each axis.*/
    void CreateCells(int ngrid)
    {
      if (fNgrid) return;
      fNgrid  = ngrid;
      fNcells = ngrid * ngrid * ngrid;
      fGrid   = new GridCell *[fNcells];
      for (int i = 0; i < fNcells; ++i)
        fGrid[i] = new GridCell();
    }

    /// Retrieve a cell by its index triplet on (x,y,z).
    /** @param ind array of cell indices on (x,y,z)
        @return GridCell corresponding to the given indices
    */
    VECCORE_ATT_HOST_DEVICE
    VECGEOM_FORCE_INLINE
    GridCell *GetCell(int ind[3]) { return fGrid[fNgrid * fNgrid * ind[0] + fNgrid * ind[1] + ind[2]]; }

    /// Retrieve a cell containing a space point and fill its indices triplet.
    /** @param[in]  point Space point
        @param[out] ind   Triplet of indices of the cell containing the point
        @return           Grid cell containing the space point
    */
    VECCORE_ATT_HOST_DEVICE
    VECGEOM_FORCE_INLINE
    GridCell *GetCell(Vector3D<T> const &point, int ind[3])
    {
      Vector3D<T> ratios = (point - fMinExtent) * fInvExtSize;
      for (int i = 0; i < 3; ++i) {
        ind[i] = ratios[i] * fNgrid;
        ind[i] = Max(ind[i], 0);
        ind[i] = Min(ind[i], fNgrid - 1);
      }
      return (GetCell(ind));
    }
  };

public:
  bool fSolidClosed   = false;   ///< Closure of the solid
  T fCubicVolume      = 0;       ///< Cubic volume
  T fSurfaceArea      = 0;       ///< Surface area
  GridHelper *fHelper = nullptr; ///< Grid helper
  Vector3D<T> fMinExtent;        ///< Minimum extent
  Vector3D<T> fMaxExtent;        ///< Maximum extent
  Vector3D<T> fInvExtSize;       ///< Inverse extent size
  Vector3D<T> fTestDir;          ///< Test direction for Inside function

  FacetBVH_t *fBVH = nullptr; ///< The BVH structure

  vector_t<Vector3D<T>> fVertices; ///< Vector of unique vertices
  vector_t<Facet_t *> fFacets;     ///< Vector of triangular facets

private:
  /// Creates the navigation acceleration structure over the builder facets.
  void CreateABBoxes()
  {
    using BoxCorner_t       = ABBoxManager<Precision>::ABBox_s;
    int nfacets             = fFacets.size();
    BoxCorner_t *boxcorners = new BoxCorner_t[2 * nfacets];
    for (int i = 0; i < nfacets; ++i) {
      Vector3D<T> minExtent(InfinityLength<T>());
      Vector3D<T> maxExtent(-InfinityLength<T>());
      auto facet = fFacets[i];
      for (size_t j = 0; j < NVERT; ++j) {
        auto &vertex = facet->fVertices[j];
        minExtent    = Min(minExtent, vertex);
        maxExtent    = Max(maxExtent, vertex);
      }
      boxcorners[2 * i]     = minExtent;
      boxcorners[2 * i + 1] = maxExtent;
    }
    fBVH = new FacetBVH_t(0, boxcorners, nfacets);

    // delete boxcorners;
  }

public:
  /// Default constructor.
  VECCORE_ATT_HOST_DEVICE
  TessellatedStruct()
  {
    fMinExtent = InfinityLength<T>();
    fMaxExtent = -InfinityLength<T>();
    fHelper    = new GridHelper();
  }

  /// Destructor.
  VECCORE_ATT_HOST_DEVICE
  ~TessellatedStruct() { delete fHelper; }

  /// Adds a pre-defined facet and re-computes the extent.
  /** The vertices are added to the list of all vertices (including duplications) and the extent
    is re-adjusted.
    @param facet Pre-computed facet to be added
  */
  VECCORE_ATT_HOST_DEVICE
  void AddFacet(Facet_t *facet)
  {

    fFacets.push_back(facet);
    int ind = fHelper->fAllVert.size();
    // Add the vertices
    for (size_t i = 0; i < NVERT; ++i) {
      fHelper->fAllVert.push_back(facet->fVertices[i]);
      facet->fIndices[i] = ind + i;
      fMinExtent[0]      = Min(fMinExtent[0], facet->fVertices[i].x());
      fMinExtent[1]      = Min(fMinExtent[1], facet->fVertices[i].y());
      fMinExtent[2]      = Min(fMinExtent[2], facet->fVertices[i].z());
      fMaxExtent[0]      = Max(fMaxExtent[0], facet->fVertices[i].x());
      fMaxExtent[1]      = Max(fMaxExtent[1], facet->fVertices[i].y());
      fMaxExtent[2]      = Max(fMaxExtent[2], facet->fVertices[i].z());
    }
  }

  /// Add a non-duplicated vertex to the solid.
  /** Duplications are only checked in the grid cell containing the vertex. An index to the unique vertex is
    added to the cell, while the vertex position is added to the list fVertices.
    @param  vertex Space point to be added as vertex.
    @return Unique vertex id after removing duplicates.
  */
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  int AddVertex(Vector3D<T> const &vertex)
  {
    // Get the cell in which to add the vertex
    int ind[3];
    GridCell *cell = fHelper->GetCell(vertex, ind);
    // Loop existing vertices in the cell and check for a duplicate
    constexpr Precision tolerancesq = kTolerance * kTolerance;
    for (int ivert : cell->fArray) {
      // existing vertex?
      if ((fVertices[ivert] - vertex).Mag2() < tolerancesq) return ivert;
    }
    // Push new vertex into the tessellated structure
    int ivertnew = fVertices.size();
    fVertices.push_back(vertex);
    // Update the cell with the new vertex index
    cell->fArray.push_back(ivertnew);
    return ivertnew;
  }

  /// Retrieval of the extent of the tessellated structure.
  /** @param[out] amin Box corner having minimum coordinates
      @param[out] amax Box corner having maximum coordinates
  */
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  void Extent(Vector3D<T> &amin, Vector3D<T> &amax) const
  {
    amin = fMinExtent;
    amax = fMaxExtent;
  }

  /// Method performing all tasks for initializing and closing the data structures.

  /** The following sequence of operations is executed:
    - Creation of the grid of cells
    - Removal of duplicate facet indices and update of facets with  the unique vertex indices
    - Creation of the navigation acceleration structures
  */
  void Close()
  {
    fInvExtSize = fMaxExtent - fMinExtent;
    if (fInvExtSize[0] * fInvExtSize[1] * fInvExtSize[2] < kTolerance) {
      std::cerr << "Tessellated structure is flat - not allowed\n";
      return;
    }
    fInvExtSize          = 1. / fInvExtSize;
    fHelper->fMinExtent  = fMinExtent;
    fHelper->fMaxExtent  = fMaxExtent;
    fHelper->fInvExtSize = fInvExtSize;
    // Make a grid with ~kVecSize facets per cell
    int ngrid = 1 + size_t(Pow<T>(T(fFacets.size()), 1. / 3.));
    // Stopwatch timer;
    // timer.Start();
    fHelper->CreateCells(ngrid);
    // auto time = timer.Stop();
    // std::cerr << "CreateCells: " << time << " sec\n";

    // Loop over facets and their vertices, fill list of vertices free of
    // duplications.
    // timer.Start();
    for (auto facet : fFacets) {
      for (int ivert = 0; ivert < 3; ++ivert) {
        facet->fIndices[ivert] = AddVertex(facet->fVertices[ivert]);
      }
    }
    // time = timer.Stop();
    // std::cerr << "Remove duplicates: " << time << " sec\n";

    // Runtime navigation consumes the scalar facet array plus BVH; the old
    // vector cluster grouping was only a legacy host-side acceleration path.
    CreateABBoxes();

    // Generate a random direction non-parallel to any of the surfaces
    constexpr T tolerance(1.e-8);
    while (true) {
      RandomDirection(fTestDir);
      bool degenerate = false;
      // Loop over triangles and check that dot product is not close to 0
      for (auto facet : fFacets) {
        if (Abs(facet->fNormal.Dot(fTestDir)) < tolerance) {
          degenerate = true;
          break;
        }
      }
      // Reject the direction, not just the rest of the facet scan.
      if (!degenerate) break;
    }
    fSolidClosed = true;
  }

  /// Generate and store a random direction used for Contains and Inside navigation queries.
  /** @param[out] direction Random direction generated */
  void RandomDirection(Vector3D<Precision> &direction)
  {
    Precision phi   = RNG::Instance().uniform(0., 2. * kPi);
    Precision theta = std::acos(1. - 2. * RNG::Instance().uniform(0, 1));
    direction.x()   = std::sin(theta) * std::cos(phi);
    direction.y()   = std::sin(theta) * std::sin(phi);
    direction.z()   = std::cos(theta);
  }

  /// Method for adding a new triangular facet
  /** @param vt0      First vertex
      @param vt1      Second vertex
      @param vt2      Third vertex
      @param absolute If true then vt0, vt1 and vt2 are the vertices to be added in
        anti-clockwise order looking from the outsider. If false the vertices are relative
        to the first: vt0, vt0+vt1, vt0+vt2, in anti-clockwise order when looking from the outsider.
  */
  VECCORE_ATT_HOST_DEVICE
  bool AddTriangularFacet(Vector3D<T> const &vt0, Vector3D<T> const &vt1, Vector3D<T> const &vt2, bool absolute = true)
  {
    VECGEOM_ASSERT(NVERT == 3);
    Facet_t *facet = new Facet_t;
    bool added     = false;
    if (absolute)
      added = facet->SetVertices(vt0, vt1, vt2);
    else
      added = facet->SetVertices(vt0, vt0 + vt1, vt0 + vt1 + vt2);
    if (!added) {
      delete facet;
      return false;
    }
    AddFacet(facet);
    return true;
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
    // We should check the quadrilateral convexity to correctly define the
    // triangle facets
    // CheckConvexity()vt0, vt1, vt2, vt3, absolute);
    VECGEOM_ASSERT(NVERT <= 4);
    Facet_t *facet = new Facet_t;
    if (NVERT == 3) {
      if (absolute) {
        if (!facet->SetVertices(vt0, vt1, vt2)) {
          delete facet;
          return false;
        }
      } else {
        if (!facet->SetVertices(vt0, vt0 + vt1, vt0 + vt1 + vt2)) {
          delete facet;
          return false;
        }
      }
      AddFacet(facet);
      facet = new Facet_t;
      if (absolute) {
        if (!facet->SetVertices(vt0, vt2, vt3)) {
          delete facet;
          return false;
        }
      } else {
        if (!facet->SetVertices(vt0, vt0 + vt2, vt0 + vt2 + vt3)) {
          delete facet;
          return false;
        }
      }
      AddFacet(facet);
      return true;
    } else if (NVERT == 4) {
      // Add a single facet
      if (absolute) {
        if (!facet->SetVertices(vt0, vt1, vt2, vt3)) {
          delete facet;
          return false;
        }
      } else {
        if (!facet->SetVertices(vt0, vt0 + vt1, vt0 + vt1 + vt2, vt0 + vt1 + vt2 + vt3)) {
          delete facet;
          return false;
        }
      }
      AddFacet(facet);
      return true;
    }
    return false;
  }
}; // end class

// Lighter weight structure holding minimal data for tessellated navigation queries
// TODO: at some point this could completely replace TessellatedStruct
template <typename T = Precision>
struct TessellatedRuntimeStruct {
  constexpr static int N = 64;

private:
  /// @brief Conservatively test one facet plane against the mesh BVH.
  /// @details Nodes whose AABB lies behind the plane are pruned. The bounded
  /// traversal avoids quadratic construction cost; an unresolved query is
  /// conservatively classified as non-supporting.
  static bool SupportsSolidPlane(TessellatedStruct<3, T> const &tsl, Vector3D<T> const &anchor,
                                 Vector3D<T> const &normal)
  {
    constexpr int kMaxNodeVisits = 2048;
    unsigned int stack[FacetBVH_t::kTraversalStackSize];
    int stackSize       = 1;
    int nodesVisited    = 0;
    int facetsRemaining = 4 * kMaxNodeVisits; // Also bound oversized leaves retained at maximum BVH depth.
    stack[0]            = 0;
    auto const *nodes   = tsl.fBVH->GetNodes();
    auto const *primIds = tsl.fBVH->GetPrimId();

    while (stackSize > 0) {
      if (++nodesVisited > kMaxNodeVisits) return false;
      auto const &node = nodes[stack[--stackSize]];
      const T x        = static_cast<T>(node.bounds[normal.x() < 0 ? 0 : 1]);
      const T y        = static_cast<T>(node.bounds[normal.y() < 0 ? 2 : 3]);
      const T z        = static_cast<T>(node.bounds[normal.z() < 0 ? 4 : 5]);
      if ((Vector3D<T>(x, y, z) - anchor).Dot(normal) <= kToleranceDist<T>) continue;

      if (node.IsLeaf()) {
        if (node.count > facetsRemaining) return false;
        facetsRemaining -= node.count;
        for (int i = 0; i < node.count; ++i) {
          auto const &facet = *tsl.fFacets[primIds[node.first + i]];
          for (int vertex = 0; vertex < 3; ++vertex)
            if ((facet.fVertices[vertex] - anchor).Dot(normal) > kToleranceDist<T>) return false;
        }
      } else {
        stack[stackSize++] = node.first;
        stack[stackSize++] = node.first + 1;
      }
    }
    return true;
  }

public:
  // pointers are used here for easier copy to GPU and allocation from existing buffers
  // TODO: treat compact memory layout
  VECCORE_ATT_HOST_DEVICE
  TessellatedRuntimeStruct(size_t ntriangles, TriangularTile<T> *t, FacetBVH_t *b,
                           unsigned char *facetSupports = nullptr)
      : fTestDir(1., 1., 1.), fNFacets(ntriangles), fFacets(t), fFacetSupports(facetSupports), fBVH(b)
  {
  }

  VECCORE_ATT_HOST_DEVICE
  TessellatedRuntimeStruct() : fFacets(nullptr), fFacetSupports(nullptr), fBVH(nullptr) {}

  Vector3D<Precision> fTestDir;
  Vector3D<Precision> fMinExtent;
  Vector3D<Precision> fMaxExtent;
  Vector3D<Precision> fExtremePoints[6]; // these points are extreme points actually on the surface of the tessellated;
                                         // Can be used for quick safety upper limit calculations

  float fTestPoints_x[N]; // safety test points in SOA mode --> target vectorized use
  float fTestPoints_y[N]; //
  float fTestPoints_z[N];

  size_t fNFacets{0}; // number of facets;

  // pointer data follows
  TriangularTile<T> *fFacets;    //
  unsigned char *fFacetSupports; ///< Per-facet supporting-plane classification.
  FacetBVH_t *fBVH;              // the BVH acceleration structure

  /// @brief Return whether the outward plane of a facet supports the complete mesh.
  /// @param facetIndex Index in the compact facet array.
  /// @return False if support was not cached or could not be certified.
  VECCORE_ATT_HOST_DEVICE
  bool FacetSupports(size_t facetIndex) const { return fFacetSupports && fFacetSupports[facetIndex] != 0; }

  /// @brief Fill compact navigation data from a tessellated builder.
  /// @param tsl Closed builder supplying the facets and acceleration structure.
  /// @param cacheFacetSupports Certify globally supporting planes with bounded
  /// BVH traversal; unresolved planes are conservatively marked non-supporting.
  void InitFrom(TessellatedStruct<3, T> const &tsl, bool cacheFacetSupports = false)
  {
    fFacets        = new TriangularTile<T>[tsl.fFacets.size()];
    fNFacets       = tsl.fFacets.size();
    fFacetSupports = cacheFacetSupports ? new unsigned char[fNFacets] : nullptr;

    double minX = InfinityLength<double>();
    double maxX = -InfinityLength<double>();
    double minY = InfinityLength<double>();
    double maxY = -InfinityLength<double>();
    double minZ = InfinityLength<double>();
    double maxZ = -InfinityLength<double>();
    int minX_i  = -1;
    int maxX_i  = -1;
    int minY_i  = -1;
    int maxY_i  = -1;
    int minZ_i  = -1;
    int maxZ_i  = -1;

    for (size_t i = 0; i < tsl.fFacets.size(); ++i) {
      const auto &facet = *(tsl.fFacets[i]);
      TriangularTile<T> t(facet.fVertices[0], facet.fVertices[1], facet.fVertices[2]);
      fFacets[i] = t;
      for (int k = 0; k < 3; ++k) {
        if (facet.fVertices[k][0] < minX) {
          minX_i = i;
          minX   = facet.fVertices[0][0];
        }
        if (facet.fVertices[k][0] > maxX) {
          maxX_i = i;
          maxX   = facet.fVertices[k][0];
        }
        if (facet.fVertices[k][1] < minY) {
          minY_i = i;
          minY   = facet.fVertices[k][1];
        }
        if (facet.fVertices[k][1] > maxY) {
          maxY_i = i;
          maxY   = facet.fVertices[k][1];
        }
        if (facet.fVertices[k][2] < minZ) {
          minZ_i = i;
          minZ   = facet.fVertices[k][2];
        }
        if (facet.fVertices[k][2] > maxZ) {
          maxZ_i = i;
          maxZ   = facet.fVertices[k][2];
        }
      }
    }
    fExtremePoints[0] = (*(tsl.fFacets[minX_i])).fVertices[0];
    fExtremePoints[1] = (*(tsl.fFacets[maxX_i])).fVertices[0];
    fExtremePoints[2] = (*(tsl.fFacets[minY_i])).fVertices[0];
    fExtremePoints[3] = (*(tsl.fFacets[maxY_i])).fVertices[0];
    fExtremePoints[4] = (*(tsl.fFacets[minZ_i])).fVertices[0];
    fExtremePoints[5] = (*(tsl.fFacets[maxZ_i])).fVertices[0];

    if (cacheFacetSupports) {
      // Local fold tolerances can accumulate across the mesh: certify each
      // plane globally, never infer support from approximate local convexity.
      Vector3D<T> supportExtremePoints[6];
      for (int i = 0; i < 6; ++i)
        supportExtremePoints[i] = tsl.fVertices[0];
      for (auto const &vertex : tsl.fVertices) {
        if (vertex.x() < supportExtremePoints[0].x()) supportExtremePoints[0] = vertex;
        if (vertex.x() > supportExtremePoints[1].x()) supportExtremePoints[1] = vertex;
        if (vertex.y() < supportExtremePoints[2].y()) supportExtremePoints[2] = vertex;
        if (vertex.y() > supportExtremePoints[3].y()) supportExtremePoints[3] = vertex;
        if (vertex.z() < supportExtremePoints[4].z()) supportExtremePoints[4] = vertex;
        if (vertex.z() > supportExtremePoints[5].z()) supportExtremePoints[5] = vertex;
      }
      for (size_t i = 0; i < fNFacets; ++i) {
        auto const &facet = fFacets[i];
        bool supports     = true;
        for (int extreme = 0; supports && extreme < 6; ++extreme)
          supports = (supportExtremePoints[extreme] - facet.fVertices[0]).Dot(facet.fNormal) <= kToleranceDist<T>;
        if (supports) supports = SupportsSolidPlane(tsl, facet.fVertices[0], facet.fNormal);
        fFacetSupports[i] = static_cast<unsigned char>(supports);
      }
    }

    // fill N random surface points --> can be used for quick upper limit
    // estimation for safety
    for (int k = 0; k < TessellatedRuntimeStruct<float>::N; ++k) {
      //      auto tp = SamplePointOnSurface(); // random index
      // for now just a random triangle --> replace by true surface point
      auto facet_index = (int)(RNG::Instance().uniform(0.0, 1.0) * fNFacets);
      auto &tp         = (*(tsl.fFacets[facet_index])).fVertices[0];
      fTestPoints_x[k] = static_cast<float>(tp[0]);
      fTestPoints_y[k] = static_cast<float>(tp[1]);
      fTestPoints_z[k] = static_cast<float>(tp[2]);
    }
    fTestDir   = tsl.fTestDir;
    fMinExtent = tsl.fMinExtent;
    fMaxExtent = tsl.fMaxExtent;

    fBVH = tsl.fBVH;
  }
};

} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif
