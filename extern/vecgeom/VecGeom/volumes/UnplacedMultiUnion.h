#ifndef VECGEOM_VOLUMES_UNPLACEDMULTIUNION_H_
#define VECGEOM_VOLUMES_UNPLACEDMULTIUNION_H_
/**
 @brief Class representing an unplaced union of multiple placed solids, possibly overlapping.
 Implemented based on the class G4MultipleUnion
 @author mihaela.gheata@cern.ch
*/

#include "VecGeom/base/Cuda.h"
#include "VecGeom/base/Global.h"
#include "VecGeom/base/AlignedBase.h"
#include "VecGeom/base/Vector3D.h"
#include "VecGeom/volumes/UnplacedVolume.h"
#include "MultiUnionStruct.h"
#include "VecGeom/volumes/kernel/MultiUnionImplementation.h"
#include "VecGeom/volumes/UnplacedVolumeImplHelper.h"

namespace vecgeom {

VECGEOM_DEVICE_DECLARE_CONV(class, UnplacedMultiUnion);
VECGEOM_DEVICE_FORWARD_DECLARE(class UnplacedMultiUnion;);

inline namespace VECGEOM_IMPL_NAMESPACE {

class UnplacedMultiUnion : public UnplacedVolumeImplHelper<MultiUnionImplementation>, public AlignedBase {

protected:
  mutable MultiUnionStruct fMultiUnion; ///< Structure storing multi-union parameters

public:
  // the constructor
  VECCORE_ATT_HOST_DEVICE
  UnplacedMultiUnion() : fMultiUnion()
  {
    fGlobalConvexity = false;
    ComputeBBox();
  }

#ifdef VECCORE_CUDA
  /**
   * Construct a non-owning runtime view over device-resident MultiUnion data.
   * @param nvolumes Number of component volumes.
   * @param volumes Device array of component placed-volume pointers.
   * @param bvh Device-resident component hierarchy.
   * @param neighbourOffsets CSR neighbour offsets.
   * @param neighbours CSR neighbour indices, or null when no components overlap.
   * @param minX Minimum X extent.
   * @param minY Minimum Y extent.
   * @param minZ Minimum Z extent.
   * @param maxX Maximum X extent.
   * @param maxY Maximum Y extent.
   * @param maxZ Maximum Z extent.
   */
  VECCORE_ATT_DEVICE
  UnplacedMultiUnion(size_t nvolumes, VPlacedVolume const **volumes, MultiUnionStruct::BVHStructure *bvh,
                     size_t *neighbourOffsets, size_t *neighbours, Precision minX, Precision minY, Precision minZ,
                     Precision maxX, Precision maxY, Precision maxZ)
      : fMultiUnion(nvolumes, volumes, bvh, neighbourOffsets, neighbours, minX, minY, minZ, maxX, maxY, maxZ)
  {
    fGlobalConvexity = false;
  }
#endif

  VECCORE_ATT_HOST
  VECGEOM_FORCE_INLINE
  void AddNode(VUnplacedVolume const *volume, Transformation3D const &transform)
  {
    LogicalVolume *lvol = new LogicalVolume(volume);
    VPlacedVolume *pvol = lvol->Place(new Transformation3D(transform));
    fMultiUnion.AddNode(pvol);
  }

  VECCORE_ATT_HOST
  VECGEOM_FORCE_INLINE
  void AddNode(VPlacedVolume const *volume) { fMultiUnion.AddNode(volume); }

  VECCORE_ATT_HOST
  void Close()
  {
    fMultiUnion.Close();
    ComputeBBox();
  }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  VPlacedVolume const *GetNode(size_t i) const { return fMultiUnion.fVolumes[i]; }

  /**
   * Replace one component pointer after placed-volume storage is compacted.
   * @param i Component index.
   * @param volume Replacement placed-volume pointer.
   */
  VECCORE_ATT_HOST
  VECGEOM_FORCE_INLINE
  void SetNode(size_t i, VPlacedVolume const *volume) { fMultiUnion.fVolumes[i] = volume; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  virtual ESolidType GetType() const override { return ESolidType::multiunion; }

  VECCORE_ATT_HOST_DEVICE
  MultiUnionStruct const &GetStruct() const { return fMultiUnion; }

  VECCORE_ATT_HOST_DEVICE
  bool Normal(Vector3D<Precision> const &point, Vector3D<Precision> &normal) const override;

  Precision Capacity() const override;

  Precision SurfaceArea() const override;

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  size_t GetNumberOfSolids() const { return fMultiUnion.fVolumes.size(); }

  VECCORE_ATT_HOST_DEVICE
  void Extent(Vector3D<Precision> &aMin, Vector3D<Precision> &aMax) const override
  {
    aMin = fMultiUnion.fMinExtent;
    aMax = fMultiUnion.fMaxExtent;
  }

  Vector3D<Precision> SamplePointOnSurface() const override;

  std::string GetEntityType() const { return "MultiUnion"; }

  VECCORE_ATT_HOST_DEVICE
  virtual void Print() const override {};

  virtual void Print(std::ostream & /*os*/) const override {};

  VECCORE_ATT_DEVICE
  static VPlacedVolume *Create(LogicalVolume const *const logical_volume, Transformation3D const *const transformation,
#ifdef VECCORE_CUDA
                               const int id, const int copy_no, const int child_id,
#endif
                               VPlacedVolume *const placement = NULL);

#ifdef VECGEOM_CUDA_INTERFACE
  virtual size_t DeviceSizeOf() const override { return DevicePtr<cuda::UnplacedMultiUnion>::SizeOf(); }
  virtual DevicePtr<cuda::VUnplacedVolume> CopyToGpu() const override;
  virtual DevicePtr<cuda::VUnplacedVolume> CopyToGpu(DevicePtr<cuda::VUnplacedVolume> const gpu_ptr) const override;
#endif

private:
  VECCORE_ATT_DEVICE
  virtual VPlacedVolume *SpecializedVolume(LogicalVolume const *const volume,
                                           Transformation3D const *const transformation,
#ifdef VECCORE_CUDA
                                           const int id, const int copy_no, const int child_id,
#endif
                                           VPlacedVolume *const placement = NULL) const override;

}; // End class
} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif /* VECGEOM_VOLUMES_UNPLACEDMULTIUNION_H_ */
