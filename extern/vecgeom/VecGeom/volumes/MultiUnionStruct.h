#ifndef VECGEOM_VOLUMES_MULTIUNIONSTRUCT_H_
#define VECGEOM_VOLUMES_MULTIUNIONSTRUCT_H_

#include "VecGeom/base/Global.h"
#include "VecGeom/base/IndexedBVH.h"
#include "VecGeom/management/ABBoxManager.h"

#include <utility>
#include <vector>

namespace vecgeom {

inline namespace VECGEOM_IMPL_NAMESPACE {

/**
 @brief Struct containing placed volume components representing a multiple union solid
 @author mihaela.gheata@cern.ch
*/

struct MultiUnionStruct {
  template <typename U>
  using vector_t = vecgeom::Vector<U>;

  using BVHReal_t    = float;
  using BVHStructure = IndexedBVH<BVHReal_t>;

  vector_t<VPlacedVolume const *> fVolumes; ///< Component placed volumes
  BVHStructure *fNavHelper = nullptr;       ///< Component AABB hierarchy

  Vector3D<Precision> fMinExtent;           ///< Minimum extent
  Vector3D<Precision> fMaxExtent;           ///< Maximum extent
  mutable Precision fCapacity    = -1;      ///< Capacity of the multiple union
  mutable Precision fSurfaceArea = -1;      ///< Surface area of the multiple union
  size_t *fNeighbourOffsets      = nullptr; ///< Component-to-neighbour range offsets
  size_t *fNeighbours            = nullptr; ///< Compact overlapping-component indices

  VECCORE_ATT_HOST_DEVICE
  MultiUnionStruct()
  {
    fMinExtent.Set(kInfLength);
    fMaxExtent.Set(-kInfLength);
  }

#ifdef VECCORE_CUDA
  /**
   * Construct a non-owning runtime view over already resident device data.
   * @param nvolumes Number of component volumes.
   * @param volumes Device array of component placed-volume pointers.
   * @param bvh Device-resident component hierarchy.
   * @param neighbourOffsets CSR offsets with @p nvolumes plus one entries.
   * @param neighbours CSR component-index payload, or null when no components overlap.
   * @param minX Minimum X extent.
   * @param minY Minimum Y extent.
   * @param minZ Minimum Z extent.
   * @param maxX Maximum X extent.
   * @param maxY Maximum Y extent.
   * @param maxZ Maximum Z extent.
   */
  VECCORE_ATT_DEVICE
  MultiUnionStruct(size_t nvolumes, VPlacedVolume const **volumes, BVHStructure *bvh, size_t *neighbourOffsets,
                   size_t *neighbours, Precision minX, Precision minY, Precision minZ, Precision maxX, Precision maxY,
                   Precision maxZ)
      : fVolumes(volumes, static_cast<int>(nvolumes)), fNavHelper(bvh), fMinExtent(minX, minY, minZ),
        fMaxExtent(maxX, maxY, maxZ), fNeighbourOffsets(neighbourOffsets), fNeighbours(neighbours)
  {
  }
#else
  MultiUnionStruct(MultiUnionStruct const &)            = delete;
  MultiUnionStruct &operator=(MultiUnionStruct const &) = delete;
#endif

  VECCORE_ATT_HOST_DEVICE
  ~MultiUnionStruct()
  {
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
    delete fNavHelper;
    delete[] fNeighbourOffsets;
    delete[] fNeighbours;
#endif
  }

  VECCORE_ATT_HOST
  void AddNode(VPlacedVolume const *volume)
  {
    Vector3D<Precision> amin, amax;
    ABBoxManager<Precision>::ComputeABBox(volume, &amin, &amax);
    fMinExtent.Set(Min(fMinExtent.x(), amin.x()), Min(fMinExtent.y(), amin.y()), Min(fMinExtent.z(), amin.z()));
    fMaxExtent.Set(Max(fMaxExtent.x(), amax.x()), Max(fMaxExtent.y(), amax.y()), Max(fMaxExtent.z(), amax.z()));
    fVolumes.push_back(volume);
  }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  static bool ABBoxOverlap(Vector3D<Precision> const &amin1, Vector3D<Precision> const &amax1,
                           Vector3D<Precision> const &amin2, Vector3D<Precision> const &amax2)
  {
    // Check if two aligned boxes overlap
    if ((amax1 - amin2).Min() < -kTolerance || (amax2 - amin1).Min() < -kTolerance) return false;
    return true;
  }

  VECCORE_ATT_HOST
  void Close()
  {
    using BoxCorner_t   = ABBoxManager<Precision>::ABBox_s;
    const size_t nboxes = fVolumes.size();
    VECGEOM_VALIDATE(nboxes > 0, << "Cannot close a MultiUnion with no component solids!");

    std::vector<BoxCorner_t> boxcorners(2 * nboxes);
    for (size_t i = 0; i < nboxes; ++i)
      ABBoxManager<Precision>::ComputeABBox(fVolumes[i], &boxcorners[2 * i], &boxcorners[2 * i + 1]);

    delete fNavHelper;
    fNavHelper = new BVHStructure(0, boxcorners.data(), static_cast<int>(nboxes));

    // Preserve the exact legacy overlap graph in a device-copyable CSR layout.
    std::vector<std::pair<size_t, size_t>> edges;
    std::vector<size_t> counts(nboxes, 0);
    for (size_t i = 0; i < nboxes - 1; ++i) {
      for (size_t j = i + 1; j < nboxes; ++j) {
        if (ABBoxOverlap(boxcorners[2 * i], boxcorners[2 * i + 1], boxcorners[2 * j], boxcorners[2 * j + 1])) {
          edges.emplace_back(i, j);
          ++counts[i];
          ++counts[j];
        }
      }
    }

    delete[] fNeighbourOffsets;
    delete[] fNeighbours;
    fNeighbourOffsets    = new size_t[nboxes + 1];
    fNeighbourOffsets[0] = 0;
    for (size_t i = 0; i < nboxes; ++i)
      fNeighbourOffsets[i + 1] = fNeighbourOffsets[i] + counts[i];

    const size_t nneighbours = fNeighbourOffsets[nboxes];
    fNeighbours              = nneighbours ? new size_t[nneighbours] : nullptr;
    std::vector<size_t> next(fNeighbourOffsets, fNeighbourOffsets + nboxes);
    for (auto const &edge : edges) {
      fNeighbours[next[edge.first]++]  = edge.second;
      fNeighbours[next[edge.second]++] = edge.first;
    }
  }

}; // End struct

} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif /* VECGEOM_VOLUMES_MULTIUNIONSTRUCT_H_ */
