/// \file BVHsurfCreator.h
/// Build and optionally dump IndexedBVH hierarchies for surface navigation.

#ifndef VECGEOM_SURFACES_BVHCREATOR_H_
#define VECGEOM_SURFACES_BVHCREATOR_H_

#include "VecGeom/base/IndexedBVH.h"
#include "VecGeom/surfaces/Model.h"

#include <cstdint>
#include <fstream>
#include <stdexcept>

namespace vgbrep {
namespace bvh {

/**
 * @brief Build a surface-navigation hierarchy from min/max corner pairs.
 *
 * The default IndexedBVH sweep-SAH policy is used because an accepted surface
 * candidate performs an exact surface query. Assignment transfers the newly
 * built arrays into the preallocated SurfData slot.
 */
template <typename Real_t>
void InitBVH(int volumeId, vecgeom::IndexedBVH<typename SurfData<Real_t>::Real_b> &bvh,
             vecgeom::Vector3D<typename SurfData<Real_t>::Real_b> const *boxes, int nBoxes)
{
  using Real_b = typename SurfData<Real_t>::Real_b;
  if (nBoxes <= 0) throw std::logic_error("Cannot construct IndexedBVH for a volume with no surfaces!");
  bvh = vecgeom::IndexedBVH<Real_b>(volumeId, boxes, nBoxes);
}

/**
 * @brief Write one IndexedBVH in a compact, versioned diagnostic format.
 *
 * The header contains `IBVH`, format version, scalar width, root id, primitive
 * count, and node count. It is followed by the primitive permutation, then
 * `(first,count,bounds[6])` for every node and six double bounds per primitive.
 */
template <typename Real_t>
void DumpBVH(vecgeom::IndexedBVH<Real_t> const &bvh, const char *filename)
{
  std::ofstream stream(filename, std::ios::binary);
  if (!stream.is_open()) return;

  constexpr std::uint32_t magic   = 0x48564249u; // "IBVH" in little-endian byte order.
  constexpr std::uint32_t version = 1;
  const std::uint32_t scalarBytes = sizeof(Real_t);
  const std::int32_t rootId       = bvh.GetRootId();
  const std::int32_t nPrimitives  = bvh.GetRootNChild();
  const std::int32_t nNodes       = bvh.GetNNodes();
  stream.write(reinterpret_cast<char const *>(&magic), sizeof(magic));
  stream.write(reinterpret_cast<char const *>(&version), sizeof(version));
  stream.write(reinterpret_cast<char const *>(&scalarBytes), sizeof(scalarBytes));
  stream.write(reinterpret_cast<char const *>(&rootId), sizeof(rootId));
  stream.write(reinterpret_cast<char const *>(&nPrimitives), sizeof(nPrimitives));
  stream.write(reinterpret_cast<char const *>(&nNodes), sizeof(nNodes));
  stream.write(reinterpret_cast<char const *>(bvh.GetPrimId()), nPrimitives * sizeof(int));

  for (int nodeId = 0; nodeId < nNodes; ++nodeId) {
    auto const &node = bvh.GetNodes()[nodeId];
    stream.write(reinterpret_cast<char const *>(&node.first), sizeof(node.first));
    stream.write(reinterpret_cast<char const *>(&node.count), sizeof(node.count));
    double bounds[6];
    for (int i = 0; i < 6; ++i)
      bounds[i] = static_cast<double>(node.bounds[i]);
    stream.write(reinterpret_cast<char const *>(bounds), sizeof(bounds));
  }

  for (int primitive = 0; primitive < nPrimitives; ++primitive) {
    auto const &box  = bvh.GetAABBs()[primitive];
    double bounds[6] = {static_cast<double>(box.Min()[0]), static_cast<double>(box.Min()[1]),
                        static_cast<double>(box.Min()[2]), static_cast<double>(box.Max()[0]),
                        static_cast<double>(box.Max()[1]), static_cast<double>(box.Max()[2])};
    stream.write(reinterpret_cast<char const *>(bounds), sizeof(bounds));
  }
}

} // namespace bvh
} // namespace vgbrep

#endif
