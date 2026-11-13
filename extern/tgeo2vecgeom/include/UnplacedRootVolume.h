/// @file unplaced_root_volume.h
/// @author Johannes de Fine Licht (johannes.definelicht@cern.ch)

#pragma once

#include "VecGeomCompat.h"

#include "VecGeom/base/Global.h"
#include "VecGeom/volumes/UnplacedVolume.h"
#include <TGeoBBox.h>
#include <TGeoShape.h>

namespace vecgeom {
namespace tgeo2vecgeom {

class UnplacedRootVolume : public VUnplacedVolume {

private:
  UnplacedRootVolume(const UnplacedRootVolume &);            // Not implemented
  UnplacedRootVolume &operator=(const UnplacedRootVolume &); // Not implemented

  TGeoShape const *fRootShape;

public:
  UnplacedRootVolume(TGeoShape const *const rootShape) : fRootShape(rootShape){};

  virtual ~UnplacedRootVolume() = default;

  TGeoShape const *GetRootShape() const { return fRootShape; }

  /// Bounding box of the underlying TGeo shape. Every TGeoShape derives from
  /// TGeoBBox, so this works for any shape we may be wrapping. Mandatory from
  /// VecGeom 2 on (pure virtual there); in VecGeom 1 the base implementation
  /// throws, so overriding it is an improvement rather than a change.
  void Extent(Vector3D<Precision> &aMin, Vector3D<Precision> &aMax) const override
  {
    auto const *box    = static_cast<TGeoBBox const *>(fRootShape);
    auto const *origin = box->GetOrigin();
    aMin.Set(origin[0] - box->GetDX(), origin[1] - box->GetDY(), origin[2] - box->GetDZ());
    aMax.Set(origin[0] + box->GetDX(), origin[1] + box->GetDY(), origin[2] + box->GetDZ());
  }

  /// Not available from TGeoShape in a form we could forward to, so report
  /// "no normal computed". Pure virtual from VecGeom 2 on.
  bool Normal(Vector3D<Precision> const & /*pos*/, Vector3D<Precision> & /*normal*/) const override { return false; }

  bool Contains(Vector3D<Precision> const &p) const override { return fRootShape->Contains(&Vector3D<double>(p)[0]); }

  EnumInside Inside(Vector3D<Precision> const &point) const override
  {
    return Contains(point) ? static_cast<EnumInside>(EInside::kInside) : static_cast<EnumInside>(EInside::kOutside);
  }

  Precision DistanceToIn(Vector3D<Precision> const &position, Vector3D<Precision> const &direction,
                         const Precision /*stepMax*/) const override
  {
    return GetRootShape()->DistFromOutside(&Vector3D<double>(position)[0], &Vector3D<double>(direction)[0], 3);
  }

  Precision DistanceToOut(Vector3D<Precision> const &position, Vector3D<Precision> const &direction,
                          const Precision /*stepMax*/) const override
  {
    return GetRootShape()->DistFromInside(&Vector3D<double>(position)[0], &Vector3D<double>(direction)[0], 3);
  }

  // NOTE: VecGeom 2 adds SurfaceHitView-taking overloads of the two functions
  // above. They are not pure, and their default bodies clear the hit info and
  // delegate to the three-argument form we do override, which is exactly what a
  // TGeo-backed shape can offer -- so there is nothing to add here.

  Precision SafetyToOut(Vector3D<Precision> const &position) const override
  {
    return GetRootShape()->Safety(&Vector3D<double>(position)[0], true);
  }

  Precision SafetyToIn(Vector3D<Precision> const &position) const override
  {
    return GetRootShape()->Safety(&Vector3D<double>(position)[0], false);
  }

  Precision Capacity() const override { return GetRootShape()->Capacity(); }

  Precision SurfaceArea() const override { return 0.; /*GetRootShape()->SurfaceArea();*/ }

  int MemorySize() const override { return sizeof(*this); }

  void Print() const override;

  void Print(std::ostream &os) const override;

private:
  virtual VPlacedVolume *SpecializedVolume(LogicalVolume const *const volume,
                                           Transformation3D const *const transformation,
#if !TGEO2VECGEOM_VECGEOM_V2
                                           const TranslationCode trans_code, const RotationCode rot_code,
#endif
                                           VPlacedVolume *const placement = NULL) const override;
};

} // End namespace tgeo2vecgeom
} // End namespace vecgeom
