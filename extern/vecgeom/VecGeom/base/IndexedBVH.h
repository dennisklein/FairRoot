/// @file IndexedBVH.h
/// @brief Compact host-built, host/device-queryable BVH with explicit child indices.
///
/// Started June 2026 by Sandro Wenzel (sandro.wenzel@cern.ch) and Tristan Wenzel
/// (tristan.wenzel@cern.ch).

#ifndef VECGEOM_BASE_INDEXEDBVH_H_
#define VECGEOM_BASE_INDEXEDBVH_H_

#include "VecGeom/base/AABB.h"
#include "VecGeom/base/Assert.h"
#include "VecGeom/base/Config.h"
#include "VecGeom/base/Cuda.h"

#include <algorithm>
#include <cstddef>
#include <cmath>
#include <cstring>
#include <limits>
#include <new>
#include <numeric>
#include <type_traits>
#include <utility>
#include <vector>

namespace vecgeom {
namespace cuda {
template <typename Real_t>
class IndexedBVH;
template <typename Real_t>
struct IndexedBVHNode;
} // namespace cuda
VECGEOM_DEVICE_DECLARE_CONV_TEMPLATE(class, IndexedBVH, typename);
VECGEOM_DEVICE_DECLARE_CONV_TEMPLATE(struct, IndexedBVHNode, typename);
inline namespace VECGEOM_IMPL_NAMESPACE {

/**
 * @brief One node in the flat explicit-index hierarchy.
 *
 * This type remains top-level so host code can name
 * `cuda::IndexedBVHNode<Real_t>` while `cuda::IndexedBVH` is forward-declared.
 *
 * A leaf has `count > 0`; `first` indexes the primitive permutation. An
 * internal node has `count == 0`; its children are `first` and `first + 1`.
 */
template <typename Real_t>
struct IndexedBVHNode {
  Real_t bounds[6]; ///< `[min_x,max_x,min_y,max_y,min_z,max_z]`, optimized for slab endpoint loads.
  int first{0};     ///< Primitive-permutation offset for a leaf, left-child node index otherwise.
  int count{0};     ///< Number of leaf primitives; zero identifies an internal node.

  /** Return whether this node is a leaf. */
  VECCORE_ATT_HOST_DEVICE
  bool IsLeaf() const { return count > 0; }

  /** Reconstruct this node's AABB from its interleaved bounds. */
  VECCORE_ATT_HOST_DEVICE
  AABB<Real_t> Bounds() const
  {
    return AABB<Real_t>(Vector3D<Real_t>(bounds[0], bounds[2], bounds[4]),
                        Vector3D<Real_t>(bounds[1], bounds[3], bounds[5]));
  }

  /** Store @p box with a one-ULP outward bias in the hierarchy precision. */
  VECCORE_ATT_HOST
  void SetBounds(const AABB<Real_t> &b)
  {
    const auto lo  = b.Min();
    const auto hi  = b.Max();
    const auto inf = std::numeric_limits<Real_t>::infinity();
    // One stored ULP covers query-point narrowing without adding work in squared-Safety traversal.
    bounds[0] = std::nextafter(lo[0], -inf);
    bounds[1] = std::nextafter(hi[0], inf);
    bounds[2] = std::nextafter(lo[1], -inf);
    bounds[3] = std::nextafter(hi[1], inf);
    bounds[4] = std::nextafter(lo[2], -inf);
    bounds[5] = std::nextafter(hi[2], inf);
  }
};

/**
 * @brief Compact explicit-index BVH for primitive or mapped-item acceleration.
 *
 * @details A binary tree is stored in one compact node array. For example:
 *
 * @verbatim
 *                         node 0: inner, first=1
 *                         /                    \
 *        node 1: leaf, first=0,count=2    node 2: inner, first=3
 *                                         /                    \
 *                        node 3: leaf, first=2,count=1    node 4: leaf, first=3,count=2
 *
 * node index:       0       1       2       3       4
 * fNodes:         inner    leaf    inner    leaf    leaf
 * first:            1       0       3       2       3
 * count:            0       2       0       1       2
 *
 * fPrimId:         [4, 1, 0, 3, 2]
 *                    ^-- node 1 uses slots [0,2): box ids 4 and 1
 *
 * box id 4 -> fAABBs[4] -> optional fBoxToDaughter[4]
 * @endverbatim
 *
 * An inner node's children are `first` and `first + 1`. A leaf owns `count`
 * consecutive slots in `fPrimId`, starting at `first`. Each slot gives the box
 * id used to access `fAABBs`. `fBoxToDaughter` is optional and allows several
 * boxes to map to one item.
 * Explicit child indices represent unbalanced trees without empty implicit-heap
 * slots; allocating siblings together preserves contiguous child loads.
 *
 * The host builder tests sweep-SAH splits on all axes. It keeps a leaf when the
 * split cost is not lower, unless the leaf cap forces a longest-axis median
 * split. Query owners traverse the explicit node indices with stacks bounded by
 * `BVH_MAX_DEPTH` and use the inlined node/primitive AABB predicates below.
 *
 * CPU and CUDA use the same arrays. The CPU slab test uses fixed endpoint loads
 * and min/max; CUDA selects endpoints with direction-octant bits.
 *
 * @tparam Real_t Scalar type stored in nodes and AABBs. Bounds are rounded
 *                    outward when converted to this type.
 */
template <typename Real_t>
class IndexedBVH {
public:
  /** @name Runtime hierarchy access
   *  These host/device declarations expose the compact arrays, traversal limits, and optional box-to-item mapping.
   *  @{ */

  /** Maximum construction depth; also bounds every fixed traversal stack. */
  static constexpr int BVH_MAX_DEPTH = 32;
  /// A depth-N binary DFS may hold N deferred siblings plus the next node.
  static constexpr int kTraversalStackSize = BVH_MAX_DEPTH + 1;

  /// Node representation stored in the flat hierarchy array.
  using Node = IndexedBVHNode<Real_t>;

  IndexedBVH() = default;

  /// Opaque owner id supplied at construction.
  VECCORE_ATT_HOST_DEVICE int GetRootId() const { return fRootId; }
  /// Number of primitive boxes in the primary hierarchy, including split boxes.
  VECCORE_ATT_HOST_DEVICE int GetRootNChild() const { return fNPrim; }
  /// Number of mapped items, or zero for an identity-mapped hierarchy.
  VECCORE_ATT_HOST_DEVICE int GetRootNDaughter() const { return fNDaughter; }
  /// Number of live nodes in the primary hierarchy.
  VECCORE_ATT_HOST_DEVICE int GetNNodes() const { return fNNodes; }
  /// Primary flat node array.
  VECCORE_ATT_HOST_DEVICE const Node *GetNodes() const { return fNodes; }
  /// Primary leaf-slot to primitive-box permutation.
  VECCORE_ATT_HOST_DEVICE const int *GetPrimId() const { return fPrimId; }
  /// Primary primitive AABBs indexed by primitive-box id.
  VECCORE_ATT_HOST_DEVICE const AABB<Real_t> *GetAABBs() const { return fAABBs; }
  /// Optional primitive-box to mapped-item table.
  VECCORE_ATT_HOST_DEVICE const int *GetBoxToDaughter() const { return fBoxToDaughter; }
  /// Whether primitive boxes are mapped to a smaller item set.
  VECCORE_ATT_HOST_DEVICE bool HasDaughterMap() const { return fBoxToDaughter != nullptr; }
  /// Resolve a primitive-box id to its mapped item id, or return it unchanged.
  VECCORE_ATT_HOST_DEVICE int GetDaughterId(int boxId) const { return fBoxToDaughter ? fBoxToDaughter[boxId] : boxId; }

  /** @} */

  /** @name Runtime AABB predicates
   *  Query owners compose these inlined host/device tests into their own traversal and exact primitive reduction.
   *  @{ */

  /**
   * @brief Test one hierarchy node against a bounded ray.
   * @param nodeId Node-array index.
   * @param point Ray origin in hierarchy precision.
   * @param inverseDirection Component-wise reciprocal direction.
   * @param octant Direction-sign bits used for CUDA endpoint selection.
   * @param step Maximum accepted ray parameter.
   * @param[out] approach Non-negative entry parameter.
   * @return Whether the ray intersects the node within @p step.
   */
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE bool IntersectNode(unsigned int nodeId, Vector3D<Real_t> const &point,
                                                                  Vector3D<Real_t> const &inverseDirection,
                                                                  unsigned int octant, Real_t step,
                                                                  Real_t &approach) const
  {
    return IntersectOctant(fNodes[nodeId].bounds, point, inverseDirection, octant, step, approach);
  }

  /**
   * @brief Test one primitive AABB against a bounded ray.
   * @param boxId Primitive-box index.
   * @param point Ray origin in hierarchy precision.
   * @param inverseDirection Component-wise reciprocal direction.
   * @param step Maximum accepted ray parameter.
   * @param[out] approach Non-negative entry parameter.
   * @return Whether the ray intersects the box within @p step.
   */
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE bool IntersectBox(int boxId, Vector3D<Real_t> const &point,
                                                                 Vector3D<Real_t> const &inverseDirection, Real_t step,
                                                                 Real_t &approach) const
  {
    return fAABBs[boxId].IntersectInvDirApproach(point, inverseDirection, step, approach);
  }

  /**
   * @brief Return conservative squared Safety to one hierarchy node.
   * @param nodeId Node-array index.
   * @param point Query point in hierarchy precision.
   * @return Downward-biased squared distance to the node AABB.
   */
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE Real_t NodeSafetySqr(unsigned int nodeId,
                                                                    Vector3D<Real_t> const &point) const
  {
    return ConservativeSafetySqr(fNodes[nodeId], point);
  }

  /**
   * @brief Apply the final conservative bias to a reduced squared Safety.
   * @param safetySqr Squared Safety reduced in hierarchy precision.
   * @return Downward-biased value converted to geometry precision.
   */
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static Precision FinalizeSafetySqr(Real_t safetySqr)
  {
    return static_cast<Precision>(kSafetyShrink * safetySqr);
  }

  /**
   * @brief Return Safety from a point to one primitive AABB.
   * @param boxId Primitive-box index.
   * @param point Query point in geometry precision.
   * @return Signed AABB Safety margin; negative when @p point is inside.
   */
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE Precision BoxSafety(int boxId, Vector3D<Precision> const &point) const
  {
    return fAABBs[boxId].Safety(point);
  }

  /**
   * @brief Test whether one hierarchy node contains a point.
   * @param nodeId Node-array index.
   * @param point Query point in hierarchy precision.
   * @return Whether the node AABB contains @p point.
   */
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE bool NodeContains(unsigned int nodeId,
                                                                 Vector3D<Real_t> const &point) const
  {
    return fNodes[nodeId].Bounds().Contains(point);
  }

  /**
   * @brief Test whether one primitive AABB contains a point.
   * @param boxId Primitive-box index.
   * @param point Query point in hierarchy precision.
   * @return Whether the primitive AABB contains @p point.
   */
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE bool BoxContains(int boxId, Vector3D<Real_t> const &point) const
  {
    return fAABBs[boxId].Contains(point);
  }

  /** @} */

  /** @name Packed device representation
   *  These declarations describe and materialize non-owning views over arrays copied into a device arena.
   *  @{ */

#ifdef VECGEOM_ENABLE_CUDA
  /** Opaque byte-offset layout for one non-owning device view and its arrays. */
  struct DeviceImageLayout {
    int rootId{0};
    int nPrim{0};
    int nDaughter{0};
    int nNodes{0};
    std::size_t nodesOffset{0};
    std::size_t primIdOffset{0};
    std::size_t aabbsOffset{0};
    std::size_t boxToDaughterOffset{0};
  };
  static_assert(std::is_trivially_copyable<DeviceImageLayout>::value,
                "IndexedBVH device-image layout must remain byte-copyable");

  /**
   * @brief Construct an identity-mapped non-owning device view.
   * @param rootId Opaque owner id available to query clients.
   * @param nPrim Number of primitive boxes.
   * @param nNodes Number of hierarchy nodes.
   * @param nodes Node array.
   * @param primId Primitive-id permutation.
   * @param aabbs Primitive AABB array.
   */
  VECCORE_ATT_HOST_DEVICE
  IndexedBVH(int rootId, int nPrim, int nNodes, Node *nodes, int *primId, AABB<Real_t> *aabbs)
      : IndexedBVH(rootId, nPrim, 0, nNodes, nodes, primId, aabbs, nullptr)
  {
  }

  /**
   * @brief Construct a non-owning device view over already resident arrays.
   *
   * @details This constructor is called by `CopyToGpu()` through placement
   * construction. Device code must not release any supplied buffer.
   *
   * @param rootId Opaque owner id available to query clients.
   * @param nPrim Number of primitive boxes in the primary hierarchy.
   * @param nDaughter Number of mapped items, or zero for identity mapping.
   * @param nNodes Number of nodes in the primary hierarchy.
   * @param nodes Primary node array.
   * @param primId Primary primitive-id permutation.
   * @param aabbs Primary primitive AABBs.
   * @param boxToDaughter Optional primitive-box to item map.
   */
  VECCORE_ATT_HOST_DEVICE
  IndexedBVH(int rootId, int nPrim, int nDaughter, int nNodes, Node *nodes, int *primId, AABB<Real_t> *aabbs,
             int *boxToDaughter)
      : fRootId(rootId), fNPrim(nPrim), fNDaughter(nDaughter), fNNodes(nNodes), fNodeCapacity(nNodes), fNodes(nodes),
        fPrimId(primId), fAABBs(aabbs), fBoxToDaughter(boxToDaughter)
  {
  }

  /**
   * @brief Replace one packed descriptor with a non-owning view over its arena.
   * @param arena Base address containing the descriptor and copied arrays.
   * @param objectOffset Byte offset of the descriptor/object slot.
   * @param target Optional separate address for the resulting view; the descriptor slot is used when null.
   */
  VECCORE_ATT_HOST_DEVICE
  static void MaterializeDeviceImage(char *arena, std::size_t objectOffset, IndexedBVH<Real_t> *target = nullptr)
  {
    DeviceImageLayout layout;
    memcpy(&layout, arena + objectOffset, sizeof(layout));

    auto *bvh = target ? target : reinterpret_cast<IndexedBVH<Real_t> *>(arena + objectOffset);
    if (layout.nPrim == 0) {
      new (bvh) IndexedBVH<Real_t>();
      return;
    }

    auto *nodes  = reinterpret_cast<Node *>(arena + layout.nodesOffset);
    auto *primId = reinterpret_cast<int *>(arena + layout.primIdOffset);
    auto *aabbs  = reinterpret_cast<AABB<Real_t> *>(arena + layout.aabbsOffset);
    auto *boxToDaughter =
        layout.boxToDaughterOffset ? reinterpret_cast<int *>(arena + layout.boxToDaughterOffset) : nullptr;

    new (bvh) IndexedBVH<Real_t>(layout.rootId, layout.nPrim, layout.nDaughter, layout.nNodes, nodes, primId, aabbs,
                                 boxToDaughter);
  }
#endif

#ifdef VECGEOM_CUDA_INTERFACE
  /**
   * @brief Reserve this BVH's arrays in a packed arena.
   * @param[in,out] nextOffset First free byte; advanced past the reserved arrays.
   * @return Opaque layout consumed by `PackDeviceImage()` and `MaterializeDeviceImage()`.
   */
  DeviceImageLayout MakeDeviceImageLayout(std::size_t &nextOffset) const;

  /**
   * @brief Write this BVH's descriptor and arrays into a host staging arena.
   * @param arena Base address of the staging arena.
   * @param arenaSize Available bytes in @p arena.
   * @param objectOffset Byte offset reserved for this BVH object.
   * @param layout Layout returned by `MakeDeviceImageLayout()` for this BVH.
   */
  void PackDeviceImage(char *arena, std::size_t arenaSize, std::size_t objectOffset,
                       DeviceImageLayout const &layout) const;

  /**
   * @brief Copy all active arrays and placement-construct a device view.
   *
   * @param addr Device address reserved for the `cuda::IndexedBVH` object.
   * @return Device pointer referring to the object constructed at @p addr.
   *
   * All active arrays share one packed allocation. The mapped-item array is
   * present only when needed. No `LogicalVolume` lookup is required.
   */
  DevicePtr<cuda::IndexedBVH<Real_t>> CopyToGpu(void *addr) const;
#endif

  /** @} */

  /** @name Host construction and ownership
   *  This section owns input validation, sweep-SAH build policy, allocation, move, and cleanup on the host.
   *  @{ */

  /** Default policy for relatively expensive work after accepting an AABB candidate. */
  static constexpr Real_t kCostTraversal = Real_t(1.0); ///< Relative internal-node traversal cost.
  static constexpr Real_t kCostIntersect = Real_t(3.0); ///< Relative exact primitive-candidate cost.
  static constexpr int kMaxLeafSize      = 4;           ///< Default maximum non-improving leaf size.

  /** @brief Host-side sweep-SAH cost model; it is not stored in the hierarchy. */
  struct BuildConfig {
    /// Relative cost of descending through an internal node and testing its children.
    double traversalCost{double(kCostTraversal)};
    /// Relative cost of accepting and processing one primitive in a leaf.
    double primitiveCost{double(kCostIntersect)};
    /// Largest leaf retained when SAH finds no cost-improving split.
    int maxLeafSize{kMaxLeafSize};
  };

  /// Binary levels of default leaves retained by the fine-grained volume policy.
  static constexpr int kVolumeFineTreeLevels = 1;
  /// Largest primitive count using the fine-grained volume policy.
  static constexpr int kVolumeWideLeafThreshold = kMaxLeafSize * (1 << kVolumeFineTreeLevels);
  /// Scale factor converting the default leaf cap to the wide-leaf cap.
  static constexpr int kVolumeWideLeafFactor = 3;

  /**
   * @brief Return a wide-leaf policy for hierarchies with cheap primitive tests.
   *
   * Up to eight boxes use `{1,2,4}`; larger hierarchies use `{2,1,12}`.
   * @p primitiveCount includes mapped or split boxes.
   * The large-tree traversal cost models two child AABB tests and permits wider
   * leaves; small trees retain finer leaves for point location.
   *
   * @param primitiveCount Number of primitive boxes in the hierarchy.
   * @return Sweep-SAH costs and leaf cap selected for this hierarchy size.
   */
  static constexpr BuildConfig VolumeNavigationBuildConfig(int primitiveCount)
  {
    return primitiveCount <= kVolumeWideLeafThreshold ? BuildConfig{1.0, 2.0, kMaxLeafSize}
                                                      : BuildConfig{2.0, 1.0, kVolumeWideLeafFactor * kMaxLeafSize};
  }

  /**
   * @brief Construct an identity-mapped hierarchy with the default policy.
   *
   * @param rootId Id of the volume/solid this BVH was built for (opaque to this class).
   * @param ptrAABB Array of @p nChild (min,max) corner pairs, one per primitive.
   * @param nChild Number of primitives.
   */
  template <typename InputReal_t>
  VECCORE_ATT_HOST IndexedBVH(int rootId, Vector3D<InputReal_t> const *ptrAABB, int nChild)
      : IndexedBVH(rootId, ptrAABB, nChild, BuildConfig{})
  {
  }

  /**
   * @brief Construct an identity-mapped hierarchy with an explicit host policy.
   *
   * @param rootId Id of the volume/solid this BVH was built for.
   * @param ptrAABB Array of @p nChild (min,max) corner pairs.
   * @param nChild Number of primitive boxes.
   * @param config Host-side construction policy.
   */
  template <typename InputReal_t>
  VECCORE_ATT_HOST IndexedBVH(int rootId, Vector3D<InputReal_t> const *ptrAABB, int nChild, BuildConfig const &config)
      : IndexedBVH(rootId, ptrAABB, nChild, nullptr, 0, config)
  {
  }

  /**
   * @brief Construct a hierarchy whose boxes may map to fewer logical items.
   *
   * @details The optional map must be dense in `[0,nDaughter)` and every item
   * must own at least one box. Repeated item ids are permitted and represent
   * conservative split boxes. Passing zero for @p nDaughter derives the count
   * from the largest mapped id.
   *
   * @param rootId Id of the volume/solid this BVH was built for.
   * @param ptrAABB Array of @p nChild (min,max) corner pairs.
   * @param nChild Number of primitive boxes.
   * @param boxToDaughter Optional primitive-box to item map. Null means identity.
   * @param nDaughter Number of mapped items, or zero to derive it from the map.
   * @param config Host-side construction policy.
   */
  template <typename InputReal_t>
  VECCORE_ATT_HOST IndexedBVH(int rootId, Vector3D<InputReal_t> const *ptrAABB, int nChild, int const *boxToDaughter,
                              int nDaughter = 0, BuildConfig const &config = BuildConfig{})
      : fRootId(rootId), fNPrim(nChild)
  {
    VECGEOM_VALIDATE(nChild > 0, << "Cannot construct IndexedBVH for a volume with no primitives!");
    VECGEOM_VALIDATE(config.traversalCost > 0. && config.primitiveCost > 0. && config.maxLeafSize > 0,
                     << "IndexedBVH build costs and maximum leaf size must be positive!");

    if (boxToDaughter) {
      VECGEOM_VALIDATE(nDaughter <= 0 || nDaughter <= nChild,
                       << "Cannot construct IndexedBVH with more mapped items than primitive boxes!");
      int mappedDaughters = 0;
      for (int boxId = 0; boxId < nChild; ++boxId) {
        int const daughterId = boxToDaughter[boxId];
        VECGEOM_VALIDATE(daughterId >= 0, << "Cannot construct IndexedBVH with a negative mapped item id!");
        VECGEOM_VALIDATE(daughterId < nChild, << "Cannot construct IndexedBVH with a sparse mapped-item id!");
        VECGEOM_VALIDATE(nDaughter <= 0 || daughterId < nDaughter,
                         << "Cannot construct IndexedBVH with an item id outside the configured range!");
        mappedDaughters = std::max(mappedDaughters, daughterId + 1);
      }
      if (nDaughter <= 0) nDaughter = mappedDaughters;

      std::vector<bool> seen(nDaughter, false);
      for (int boxId = 0; boxId < nChild; ++boxId)
        seen[boxToDaughter[boxId]] = true;
      VECGEOM_VALIDATE(std::all_of(seen.begin(), seen.end(), [](bool value) { return value; }),
                       << "Cannot construct IndexedBVH when a mapped item has no primitive box!");

      fNDaughter     = nDaughter;
      fBoxToDaughter = new int[nChild];
      std::copy(boxToDaughter, boxToDaughter + nChild, fBoxToDaughter);
    }

    fPrimId = new int[nChild];
    std::iota(fPrimId, fPrimId + nChild, 0);

    // Per-primitive AABBs indexed by primitive id. Corners are rounded outward when narrowing the input,
    // then expanded by one stored ULP to cover nearby query-point narrowing,
    // including when the source bounds are exactly representable.
    fAABBs = new AABB<Real_t>[nChild];
    for (int i = 0; i < nChild; ++i) {
      Vector3D<Real_t> lo(RoundOutward<true>(ptrAABB[2 * i].x()), RoundOutward<true>(ptrAABB[2 * i].y()),
                          RoundOutward<true>(ptrAABB[2 * i].z()));
      Vector3D<Real_t> hi(RoundOutward<false>(ptrAABB[2 * i + 1].x()), RoundOutward<false>(ptrAABB[2 * i + 1].y()),
                          RoundOutward<false>(ptrAABB[2 * i + 1].z()));
      const auto inf = std::numeric_limits<Real_t>::infinity();
      for (int axis = 0; axis < 3; ++axis) {
        lo[axis] = std::nextafter(lo[axis], -inf);
        hi[axis] = std::nextafter(hi[axis], inf);
      }
      fAABBs[i] = AABB<Real_t>(lo, hi);
    }

    // A binary tree with at most nChild leaves has at most 2 * nChild - 1 nodes.
    fNodeCapacity = std::max(1, 2 * nChild - 1);
    fNodes        = new Node[fNodeCapacity];
    fNNodes       = 1;
    BuildNode(0, 0, fAABBs, fPrimId, fPrimId + nChild, config);

    // The recursive builder needs worst-case 2*N-1 scratch capacity. SAH normally stops at wider
    // leaves, so discard the unused tail once construction is complete. Traversal indices remain
    // unchanged because the live nodes already occupy the contiguous prefix [0, fNNodes).
    if (fNNodes < fNodeCapacity) {
      auto *compactNodes = new Node[fNNodes];
      std::copy(fNodes, fNodes + fNNodes, compactNodes);
      delete[] fNodes;
      fNodes        = compactNodes;
      fNodeCapacity = fNNodes;
    }
  }

  IndexedBVH(IndexedBVH const &)            = delete;
  IndexedBVH &operator=(IndexedBVH const &) = delete;

  /** Transfer ownership of all host arrays. */
  VECCORE_ATT_HOST
  IndexedBVH(IndexedBVH &&other) noexcept { MoveFrom(other); }

  /** Release current host arrays and take ownership from @p other. */
  VECCORE_ATT_HOST
  IndexedBVH &operator=(IndexedBVH &&other) noexcept
  {
    if (this != &other) {
      Clear();
      MoveFrom(other);
    }
    return *this;
  }

  /** @return Bytes owned by the active host arrays, excluding the object itself. */
  std::size_t GetAllocatedSize() const
  {
    std::size_t size = static_cast<std::size_t>(fNodeCapacity) * sizeof(Node) +
                       static_cast<std::size_t>(fNPrim) * (sizeof(int) + sizeof(AABB<Real_t>));
    if (fBoxToDaughter) size += static_cast<std::size_t>(fNPrim) * sizeof(int);
    return size;
  }

  /** @return Active device-array payload bytes, excluding layout, alignment, and object storage. */
  std::size_t GetDeviceAllocatedSize() const
  {
    std::size_t size = static_cast<std::size_t>(fNNodes) * sizeof(Node) +
                       static_cast<std::size_t>(fNPrim) * (sizeof(int) + sizeof(AABB<Real_t>));
    if (fBoxToDaughter) size += static_cast<std::size_t>(fNPrim) * sizeof(int);
    return size;
  }

  ~IndexedBVH() { Clear(); }

  /** Release every host-owned node, primitive, AABB, and mapping array. */
  VECCORE_ATT_HOST
  void Clear()
  {
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
    delete[] fNodes;
    fNodes = nullptr;
    delete[] fPrimId;
    fPrimId = nullptr;
    delete[] fAABBs;
    fAABBs = nullptr;
    delete[] fBoxToDaughter;
    fBoxToDaughter = nullptr;
    fRootId        = 0;
    fNPrim         = 0;
    fNDaughter     = 0;
    fNNodes        = 0;
    fNodeCapacity  = 0;
#endif
  }

  /** @} */

  /** @name Host build diagnostics
   *  These host-only estimates characterize the completed hierarchy without participating in runtime queries.
   *  @{ */

  /**
   * @brief Estimate isotropic ray work using surface-area probabilities.
   * @param[out] nodeVisits Estimated internal and leaf node visits.
   * @param[out] primitiveBoxTests Estimated primitive AABB tests.
   * @param[out] candidateCalls Estimated exact primitive operations.
   */
  VECCORE_ATT_HOST void EstimateTraversalWork(double &nodeVisits, double &primitiveBoxTests,
                                              double &candidateCalls) const
  {
    nodeVisits = primitiveBoxTests = candidateCalls = 0.;
    if (!fNodes || !fPrimId || !fAABBs) return;

    const double rootArea = static_cast<double>(fNodes[0].Bounds().SurfaceArea());
    if (!(rootArea > 0.)) return;

    unsigned int stack[kTraversalStackSize];
    int sp   = 1;
    stack[0] = 0;
    while (sp > 0) {
      const unsigned int id    = stack[--sp];
      const double probability = id == 0 ? 1. : static_cast<double>(fNodes[id].Bounds().SurfaceArea()) / rootArea;
      nodeVisits += probability;
      if (fNodes[id].IsLeaf()) {
        primitiveBoxTests += probability * fNodes[id].count;
        for (int child = 0; child < fNodes[id].count; ++child) {
          const int boxId = fPrimId[fNodes[id].first + child];
          candidateCalls += static_cast<double>(fAABBs[boxId].SurfaceArea()) / rootArea;
        }
      } else {
        stack[sp++] = fNodes[id].first;
        stack[sp++] = fNodes[id].first + 1;
      }
    }
  }

  /**
   * @brief Estimate point-Safety work for the primary hierarchy.
   * @details Node visits use boundary-biased surface-area ratios; leaf and
   * candidate work use occupied-volume ratios.
   * @param[out] nodeVisits Estimated node visits.
   * @param[out] primitiveBoxTests Estimated primitive AABB tests.
   * @param[out] candidateCalls Estimated exact item operations.
   */
  VECCORE_ATT_HOST void EstimateSafetyWork(double &nodeVisits, double &primitiveBoxTests, double &candidateCalls) const
  {
    nodeVisits = primitiveBoxTests = candidateCalls = 0.;
    if (!fNodes || !fPrimId || !fAABBs) return;

    const auto volume = [](AABB<Real_t> const &box) {
      const auto extent = box.Size();
      return static_cast<double>(extent[0]) * static_cast<double>(extent[1]) * static_cast<double>(extent[2]);
    };
    const double rootArea   = static_cast<double>(fNodes[0].Bounds().SurfaceArea());
    const double rootVolume = volume(fNodes[0].Bounds());
    if (!(rootArea > 0.) || !(rootVolume > 0.)) return;

    unsigned int stack[kTraversalStackSize];
    int sp   = 1;
    stack[0] = 0;
    while (sp > 0) {
      const unsigned int id = stack[--sp];
      nodeVisits += id == 0 ? 1. : static_cast<double>(fNodes[id].Bounds().SurfaceArea()) / rootArea;
      if (fNodes[id].IsLeaf()) {
        const double occupancy = id == 0 ? 1. : volume(fNodes[id].Bounds()) / rootVolume;
        primitiveBoxTests += occupancy * fNodes[id].count;
        for (int child = 0; child < fNodes[id].count; ++child) {
          const int boxId = fPrimId[fNodes[id].first + child];
          candidateCalls += volume(fAABBs[boxId]) / rootVolume;
        }
      } else {
        stack[sp++] = fNodes[id].first;
        stack[sp++] = fNodes[id].first + 1;
      }
    }
  }

  /** @} */

private:
  // Runtime storage
  // The same pointer-only representation is used by owning host objects and non-owning device views.
  int fRootId{0};                ///< Opaque owner id available to query clients.
  int fNPrim{0};                 ///< Number of primitive boxes in the primary hierarchy.
  int fNDaughter{0};             ///< Number of mapped items; zero means identity mapping.
  int fNNodes{0};                ///< Live nodes in the primary hierarchy.
  int fNodeCapacity{0};          ///< Host allocation capacity; compacted to `fNNodes` after build.
  Node *fNodes{nullptr};         ///< Flat primary hierarchy.
  int *fPrimId{nullptr};         ///< Leaf-slot to primitive-box id permutation.
  AABB<Real_t> *fAABBs{nullptr}; ///< Primary AABBs indexed by primitive-box id.
  int *fBoxToDaughter{nullptr};  ///< Optional primitive-box id to mapped-item id.

  // Runtime predicate implementation
  // These helpers implement the public host/device box predicates without navigation semantics.
  /// Downward relative bias applied to squared AABB Safety in hierarchy precision.
  static constexpr Real_t kSafetyShrink = Real_t(1.) - Real_t(8.) * vecgeom::kEpsilonT<Real_t>;
  /**
   * @brief Return a conservative squared distance to one node.
   * @details Bounds are already widened. Multiplying the squared result by
   * `kSafetyShrink` also covers accumulated `Real_t` subtraction, max, square,
   * and addition rounding without a square root in the hot traversal.
   */
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static Real_t ConservativeSafetySqr(const Node &node,
                                                                                   const Vector3D<Real_t> &point)
  {
    const Real_t dx = vecgeom::Max(vecgeom::Max(node.bounds[0] - point[0], Real_t(0.)), point[0] - node.bounds[1]);
    const Real_t dy = vecgeom::Max(vecgeom::Max(node.bounds[2] - point[1], Real_t(0.)), point[1] - node.bounds[3]);
    const Real_t dz = vecgeom::Max(vecgeom::Max(node.bounds[4] - point[2], Real_t(0.)), point[2] - node.bounds[5]);
    return kSafetyShrink * (dx * dx + dy * dy + dz * dz);
  }

  /**
   * @brief Test a ray against one node AABB.
   *
   * CUDA uses @p octant to select slab endpoints. CPU loads both endpoints and
   * applies min/max. Both compute the same entry/exit interval, preserving hit
   * decisions and near-child ordering within `[0,step]`.
   *
   * @param bounds Interleaved minimum/maximum endpoints of the node AABB.
   * @param point Ray origin in the hierarchy coordinate system.
   * @param invdir Component-wise reciprocal ray direction.
   * @param octant Direction-sign bits used by the CUDA endpoint loads; ignored on the host.
   * @param step Maximum accepted ray parameter.
   * @param[out] tmin Entry distance, clamped to >= 0; used by the caller for near/far child ordering.
   * @return Whether the ray crosses @p bounds within `[0,step]`.
   */
  VECCORE_ATT_HOST_DEVICE
  static bool IntersectOctant(const Real_t (&bounds)[6], const Vector3D<Real_t> &point, const Vector3D<Real_t> &invdir,
                              unsigned int octant, Real_t step, Real_t &tmin)
  {
#ifdef VECCORE_CUDA_DEVICE_COMPILATION
    // bit_i is 0/1 so 2*axis+bit_i / 2*axis+(1-bit_i) are plain indexed loads, no select per axis.
    const unsigned int bit0 = octant & 1u;
    const unsigned int bit1 = (octant >> 1) & 1u;
    const unsigned int bit2 = (octant >> 2) & 1u;

    const Real_t t0x = (bounds[bit0] - point[0]) * invdir[0];
    const Real_t t1x = (bounds[1 - bit0] - point[0]) * invdir[0];
    const Real_t t0y = (bounds[2 + bit1] - point[1]) * invdir[1];
    const Real_t t1y = (bounds[3 - bit1] - point[1]) * invdir[1];
    const Real_t t0z = (bounds[4 + bit2] - point[2]) * invdir[2];
    const Real_t t1z = (bounds[5 - bit2] - point[2]) * invdir[2];

    using vecCore::math::Max;
    using vecCore::math::Min;
    tmin = Max(Max(Max(t0x, t0y), t0z), Real_t(0.));
    // Same outward bias on the exit t as AABB::ComputeIntersectionInvDir, for float-narrowing safety.
    const Real_t tmax = Min(Min(Min(t1x, t1y), t1z) * (Real_t(1.) + vecgeom::kToleranceDist<Real_t>), step);
    return tmin <= tmax;
#else
    // Load both endpoints from fixed addresses; min/max assigns near and far after multiplication.
    // This is algebraically identical to selecting the bounds first from the direction octant.
    (void)octant;
    const Real_t t0x = (bounds[0] - point[0]) * invdir[0];
    const Real_t t1x = (bounds[1] - point[0]) * invdir[0];
    const Real_t t0y = (bounds[2] - point[1]) * invdir[1];
    const Real_t t1y = (bounds[3] - point[1]) * invdir[1];
    const Real_t t0z = (bounds[4] - point[2]) * invdir[2];
    const Real_t t1z = (bounds[5] - point[2]) * invdir[2];

    using vecCore::math::Max;
    using vecCore::math::Min;
    tmin              = Max(Max(Max(Min(t0x, t1x), Min(t0y, t1y)), Min(t0z, t1z)), Real_t(0.));
    const Real_t tmax = Min(
        Min(Min(Max(t0x, t1x), Max(t0y, t1y)), Max(t0z, t1z)) * (Real_t(1.) + vecgeom::kToleranceDist<Real_t>), step);
    return tmin <= tmax;
#endif
  }

  // Host construction implementation
  // Ownership transfer, conservative input conversion, SAH partitioning, and recursive building remain host-only.
  VECCORE_ATT_HOST
  void MoveFrom(IndexedBVH &other) noexcept
  {
    fRootId        = std::exchange(other.fRootId, 0);
    fNPrim         = std::exchange(other.fNPrim, 0);
    fNDaughter     = std::exchange(other.fNDaughter, 0);
    fNNodes        = std::exchange(other.fNNodes, 0);
    fNodeCapacity  = std::exchange(other.fNodeCapacity, 0);
    fNodes         = std::exchange(other.fNodes, nullptr);
    fPrimId        = std::exchange(other.fPrimId, nullptr);
    fAABBs         = std::exchange(other.fAABBs, nullptr);
    fBoxToDaughter = std::exchange(other.fBoxToDaughter, nullptr);
  }

  /**
   * @brief Narrow a bound to `Real_t` without shrinking its AABB.
   *
   * @tparam isMin Round down when true and up otherwise.
   * @tparam From Source scalar type.
   * @param v Bound to convert.
   * @return The converted bound, advanced by one representable value only when
   *         the ordinary conversion rounded inward.
   */
  template <bool isMin, typename From>
  VECCORE_ATT_HOST static Real_t RoundOutward(From v)
  {
    Real_t r = static_cast<Real_t>(v);
    if (isMin ? (static_cast<From>(r) > v) : (static_cast<From>(r) < v)) {
      r = isMin ? std::nextafter(r, -std::numeric_limits<Real_t>::infinity())
                : std::nextafter(r, std::numeric_limits<Real_t>::infinity());
    }
    return r;
  }

  /**
   * @brief Lexicographically compare 3D keys starting with one axis.
   *
   * @tparam T Three-component key type.
   * @param left First key.
   * @param right Second key.
   * @param sortAxis Primary axis; the other two axes break ties cyclically.
   * @return Whether @p left precedes @p right. Identical keys remain equivalent.
   */
  template <typename T>
  VECCORE_ATT_HOST static bool less3D(const T &left, const T &right, int sortAxis)
  {
    return left[sortAxis] < right[sortAxis] ||
           (left[sortAxis] == right[sortAxis] && (left[(sortAxis + 1) % 3] < right[(sortAxis + 1) % 3] ||
                                                  (left[(sortAxis + 1) % 3] == right[(sortAxis + 1) % 3] &&
                                                   left[(sortAxis + 2) % 3] < right[(sortAxis + 2) % 3])));
  }

  /**
   * @brief Compute prefix/suffix surface areas for one sorted primitive range.
   *
   * @details The returned `[k].first` covers `[begin,begin+k)`, while
   * `[k].second` covers `[begin+k,end)`, allowing every partition cost on an
   * axis to be evaluated in linear time after sorting.
   *
   * @param boxes Primitive AABBs indexed by the values in the range.
   * @param begin First primitive-id slot.
   * @param end One-past-last primitive-id slot.
   * @return Prefix/suffix enclosing-box surface areas for all split positions.
   */
  VECCORE_ATT_HOST
  static std::vector<std::pair<double, double>> SweepSurfaceArea(AABB<Real_t> const *boxes, int const *begin,
                                                                 int const *end)
  {
    if (begin >= end) return {};

    std::vector<std::pair<double, double>> areas(std::distance(begin, end), {0., 0.});

    AABB<Real_t> box{boxes[*begin]};
    for (auto it = begin + 1; it < end; ++it) {
      areas[it - begin].first = box.SurfaceArea();
      box                     = AABB<Real_t>::Union(box, boxes[*it]);
    }

    AABB<Real_t> box2{boxes[*(end - 1)]};
    for (auto it = end; it != begin;) {
      --it;
      box2                     = AABB<Real_t>::Union(box2, boxes[*it]);
      areas[it - begin].second = box2.SurfaceArea();
    }

    return areas;
  }

  /**
   * @brief Find a cost-improving sweep-SAH partition over all three axes.
   *
   * Minimizes `c_trav + c_prim * sum(P(child) * N(child))` against the leaf
   * cost. A small extent term orders boxes with equal centroids.
   *
   * @param boxes Primitive AABBs indexed by the mutable range.
   * @param first First primitive-id slot.
   * @param last One-past-last primitive-id slot.
   * @param nodeArea Surface area of the range's enclosing AABB.
   * @param config Relative traversal/primitive costs and leaf cap.
   * @param[out] improves Whether a strictly cost-improving split was found.
   * @return First element of the second group when @p improves, otherwise @p last.
   */
  VECCORE_ATT_HOST
  static int *SweepSAHSplit(AABB<Real_t> const *boxes, int *first, int *last, Real_t nodeArea,
                            BuildConfig const &config, bool &improves)
  {
    const auto n         = std::distance(first, last);
    double bestCost      = double(n) * config.primitiveCost; // baseline: a single leaf
    int bestSplitAxis    = -1;
    int bestSplitOffset  = -1;
    const double invArea = 1.0 / double(nodeArea);

    int currentSortAxis = 0;
    auto sorter         = [boxes, &currentSortAxis](int a, int b) {
      const auto centroidA   = boxes[a].Center();
      const auto centroidB   = boxes[b].Center();
      constexpr double shift = 0.01;
      return less3D(centroidA + shift * (centroidA - boxes[a].Min()), centroidB + shift * (centroidB - boxes[b].Min()),
                    currentSortAxis);
    };

    for (int axis = 0; axis <= 2; ++axis) {
      currentSortAxis = axis;
      std::sort(first, last, sorter);

      const std::vector<std::pair<double, double>> sweep = SweepSurfaceArea(boxes, first, last);

      // k = number of elements in the left group; both groups must be non-empty.
      for (int *splitObject = first + 1; splitObject < last; ++splitObject) {
        const auto k      = std::distance(first, splitObject);
        const double cost = config.traversalCost +
                            config.primitiveCost * (sweep[k].first * invArea * k + sweep[k].second * invArea * (n - k));

        if (cost < bestCost) {
          bestCost        = cost;
          bestSplitAxis   = axis;
          bestSplitOffset = static_cast<int>(k);
        }
      }
    }

    improves = (bestSplitAxis != -1);
    if (!improves) return last;

    currentSortAxis = bestSplitAxis;
    std::sort(first, last, sorter);
    return first + bestSplitOffset;
  }

  /**
   * @brief Partition at the median centroid along the longest node axis.
   *
   * Used when the leaf cap requires a split that SAH rejects.
   *
   * @param boxes Primitive AABBs indexed by the mutable range.
   * @param first First primitive-id slot.
   * @param last One-past-last primitive-id slot.
   * @param bounds Enclosing AABB used to select the longest axis.
   * @return First primitive-id slot in the right partition.
   */
  VECCORE_ATT_HOST
  static int *MedianSplitLongestAxis(AABB<Real_t> const *boxes, int *first, int *last, AABB<Real_t> const &bounds)
  {
    const auto size = bounds.Size();
    const int axis  = (size[0] > size[2]) ? (size[0] > size[1] ? 0 : 1) : (size[1] > size[2] ? 1 : 2);
    int *mid        = first + std::distance(first, last) / 2;
    std::nth_element(first, mid, last,
                     [boxes, axis](int a, int b) { return less3D(boxes[a].Center(), boxes[b].Center(), axis); });
    return mid;
  }

  /**
   * @brief Recursively build one node for the primitive range `[first,last)`.
   * @details Stores the range union, accepts a leaf at one primitive or maximum
   * depth, otherwise applies sweep SAH. A non-improving range remains a leaf up
   * to `BuildConfig::maxLeafSize`; larger ranges use the longest-axis median.
   * Child nodes are allocated as one contiguous pair before recursion.
   *
   * @param id Node-array slot to populate.
   * @param depth Current root-relative depth.
   * @param boxes Primitive AABBs indexed by values in `[first,last)`.
   * @param first First primitive-id slot owned by this node.
   * @param last One-past-last primitive-id slot owned by this node.
   * @param config Construction cost model and maximum leaf size.
   */
  VECCORE_ATT_HOST
  void BuildNode(int id, int depth, AABB<Real_t> const *boxes, int *first, int *last, BuildConfig const &config)
  {
    Node &node  = fNodes[id];
    const int n = static_cast<int>(std::distance(first, last));

    AABB<Real_t> bbox = boxes[*first];
    for (int *it = first + 1; it != last; ++it)
      bbox = AABB<Real_t>::Union(bbox, boxes[*it]);
    node.SetBounds(bbox);

    auto makeLeaf = [&] {
      node.first = static_cast<int>(first - fPrimId);
      node.count = n;
    };

    if (n == 1 || depth >= BVH_MAX_DEPTH) return makeLeaf();

    const Real_t nodeArea = bbox.SurfaceArea();
    bool improves         = false;
    int *pivot = (nodeArea > Real_t(0.)) ? SweepSAHSplit(boxes, first, last, nodeArea, config, improves) : last;

    if (!improves) {
      if (n <= config.maxLeafSize) return makeLeaf();
      // The leaf-size cap still forces a split: fall back to a median split.
      pivot = MedianSplitLongestAxis(boxes, first, last, bbox);
      // Truly degenerate (all primitives coincide): accept an oversized leaf rather than loop forever.
      if (pivot == first || pivot == last) return makeLeaf();
    }

    const int childL = fNNodes++;
    const int childR = fNNodes++;
    node.first       = childL;
    node.count       = 0;

    BuildNode(childL, depth + 1, boxes, first, pivot, config);
    BuildNode(childR, depth + 1, boxes, pivot, last, config);
  }
};

} // namespace VECGEOM_IMPL_NAMESPACE

#ifdef VECGEOM_ENABLE_CUDA
/**
 * @brief Copy a contiguous host array of IndexedBVHs into one device arena.
 *
 * The returned pointer is both the device object-array base and the allocation
 * to release. Each object materializes its own arrays from the packed arena.
 */
template <typename Real_t>
cuda::IndexedBVH<Real_t> *CopyManyIndexedBVHToGpu(void const *hostBVHs, std::size_t count);
#endif

} // namespace vecgeom

#endif // VECGEOM_BASE_INDEXEDBVH_H_
