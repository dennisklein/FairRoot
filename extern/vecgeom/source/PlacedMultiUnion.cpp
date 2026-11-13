/// @file PlacedMultiUnion.cpp
/// @author Mihaela Gheata (mihaela.gheata@cern.ch)

#include "VecGeom/volumes/MultiUnion.h"

#ifndef VECCORE_CUDA
#ifdef VECGEOM_ROOT
#include <TGeoCompositeShape.h>
#include <TGeoBoolNode.h>
#include <TGeoMatrix.h>
#endif
#ifdef VECGEOM_GEANT4
#include <G4MultiUnion.hh>
#include <G4ThreeVector.hh>
#include <G4RotationMatrix.hh>
#include <G4Transform3D.hh>
#endif

#endif // VECCORE_CUDA

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

VECCORE_ATT_HOST_DEVICE
void PlacedMultiUnion::PrintType() const { printf("PlacedMultiUnion"); }

void PlacedMultiUnion::PrintType(std::ostream &s) const { s << "PlacedMultiUnion"; }

#ifndef VECCORE_CUDA

VPlacedVolume const *PlacedMultiUnion::ConvertToUnspecialized() const
{
  return new SimpleMultiUnion(GetLabel().c_str(), GetLogicalVolume(), GetTransformation());
}

#ifdef VECGEOM_ROOT
TGeoShape const *PlacedMultiUnion::ConvertToRoot() const
{
  auto const *multi = GetUnplacedVolume();
  const auto count  = multi->GetNumberOfSolids();
  if (count == 0) return nullptr;
  TGeoShape *result = nullptr;
  TGeoMatrix *frame = nullptr;
  for (size_t i = 0; i < count; ++i) {
    auto const *placed = multi->GetNode(i);
    auto *shape        = const_cast<TGeoShape *>(placed->ConvertToRoot());
    if (!shape) return nullptr;
    auto *transform = Transformation3D::ConvertToTGeoMatrix(*placed->GetTransformation());
    if (i == 0) {
      result = shape;
      frame  = transform;
    } else {
      result = new TGeoCompositeShape("RootMultiUnion", new TGeoUnion(result, shape, frame, transform));
      frame  = nullptr;
    }
  }
  // A single component still needs its placement expressed in the union frame.
  if (count == 1) result = new TGeoCompositeShape("RootMultiUnion", new TGeoUnion(result, result, frame, frame));
  return result;
}
#endif

#ifdef VECGEOM_GEANT4
G4VSolid const *PlacedMultiUnion::ConvertToGeant4() const
{
  G4MultiUnion *munion               = new G4MultiUnion("g4multiunion");
  const UnplacedMultiUnion *unplaced = GetUnplacedVolume();
  for (size_t i = 0; i < unplaced->GetNumberOfSolids(); ++i) {
    G4VSolid *g4solid = (G4VSolid *)unplaced->GetNode(i)->ConvertToGeant4();
    if (!g4solid) {
      delete munion;
      return nullptr;
    }
    auto trans = unplaced->GetNode(i)->GetTransformation();
    // Vector3D<double> point(1, 1, 1);
    // Vector3D<double> pnew = trans->Transform(point);
    const double *rotm = trans->Rotation();
    G4RotationMatrix g4rot(G4ThreeVector(rotm[0], rotm[3], rotm[6]), G4ThreeVector(rotm[1], rotm[4], rotm[7]),
                           G4ThreeVector(rotm[2], rotm[5], rotm[8]));
    G4Transform3D g4trans(g4rot, G4ThreeVector(trans->Translation(0), trans->Translation(1), trans->Translation(2)));
    // G4ThreeVector global(1, 1, 1);
    // G4ThreeVector local = g4trans.inverse()*G4Point3D(global);
    munion->AddNode(*g4solid, g4trans);
  }
  munion->Voxelize();
  return munion;
}
#endif

#endif // VECCORE_CUDA

} // namespace VECGEOM_IMPL_NAMESPACE

#ifdef VECCORE_CUDA

VECGEOM_DEVICE_INST_PLACED_VOLUME_ALLSPEC(SpecializedMultiUnion)

#endif

} // namespace vecgeom
