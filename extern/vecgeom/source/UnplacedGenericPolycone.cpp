/// @file UnplacedGenericPolycone.cpp
/// @author Raman Sehgal (raman.sehgal@cern.ch)

#include "VecGeom/volumes/Cone.h"
#include "VecGeom/volumes/UnplacedCoaxialCones.h"
#include "VecGeom/management/GeoManager.h"
#include "VecGeom/volumes/UnplacedGenericPolycone.h"
#include "VecGeom/management/VolumeFactory.h"
#include "VecGeom/volumes/SpecializedGenericPolycone.h"
#include "VecGeom/base/RNG.h"
#ifndef VECCORE_CUDA
#include "VecGeom/volumes/utilities/BoundingEnvelope.h"
#endif
#include <stdio.h>
#include <cmath>

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

#ifndef VECCORE_CUDA
namespace {

/// Builds a conservative union of angular frustum cells from reduced sections.
///
/// GenericPolycone sections may contain several coaxial cone pieces. Each piece
/// is already a section-local conical frustum, so the clipped extent can reuse
/// BoundingEnvelope's circumscribed angular-cell envelope without consulting
/// the original visualization contour.
///
/// @return True when at least one finite coaxial-cone envelope was added.
bool AddGenericPolyconeEnvelope(BoundingEnvelope &envelope, UnplacedGenericPolycone const &polycone)
{
  GenericPolyconeStruct<Precision> const &storage = polycone.GetStruct();
  bool added_piece                                = false;
  for (int section_index = 0; section_index < storage.GetNSections(); ++section_index) {
    GenericPolyconeSection const &section              = storage.GetSection(section_index);
    CoaxialConesStruct<Precision> const *coaxial_cones = section.fCoaxialCones;
    if (coaxial_cones == nullptr) continue;

    const Precision dz = Abs(coaxial_cones->fDz);
    if (dz <= kTolerance) continue;
    const Precision z1 = section.fShift - dz;
    const Precision z2 = section.fShift + dz;

    for (unsigned int cone_index = 0; cone_index < coaxial_cones->fConeStructVector.size(); ++cone_index) {
      ConeStruct<Precision> const *cone = coaxial_cones->fConeStructVector[cone_index];
      if (cone == nullptr) continue;
      const Precision rmax1 = cone->fRmax1;
      const Precision rmax2 = cone->fRmax2;
      envelope.AddLocalConicalFrustum(cone->fRmin1, rmax1, cone->fRmin2, rmax2, z1, z2, polycone.GetSPhi(),
                                      polycone.GetDPhi());
      added_piece = true;
    }
  }
  return added_piece;
}

} // namespace
#endif

/*
 * All the required Parametric Constructor
 */

VECCORE_ATT_HOST_DEVICE
UnplacedGenericPolycone::UnplacedGenericPolycone(Precision phiStart, // initial phi starting angle
                                                 Precision phiTotal, // total phi angle
                                                 int numRZ, // number corners in r,z space (must be an even number)
                                                 Precision const *r, // r coordinate of these corners
                                                 Precision const *z)
{
  fGlobalConvexity = false;
  fSPhi            = phiStart;
  fDPhi            = phiTotal;
  fNumRZ           = numRZ;
  for (int i = 0; i < numRZ; i++) {
    fR.push_back(r[i]);
    fZ.push_back(z[i]);
  }

  Vector<Vector2D<Precision>> rzVect;
  for (int i = 0; i < numRZ; i++) {
    rzVect.push_back(Vector2D<Precision>(r[i], z[i]));
  }
  ReducedPolycone rd(rzVect);
  Vector<Precision> zS;
  Vector<Vector<Precision>> vectOfRmin1Vect;
  Vector<Vector<Precision>> vectOfRmax1Vect;
  Vector<Vector<Precision>> vectOfRmin2Vect;
  Vector<Vector<Precision>> vectOfRmax2Vect;
  rd.GetPolyconeParameters(vectOfRmin1Vect, vectOfRmax1Vect, vectOfRmin2Vect, vectOfRmax2Vect, zS, fAMin, fAMax);

  fGenericPolycone.Set(vectOfRmin1Vect, vectOfRmax1Vect, vectOfRmin2Vect, vectOfRmax2Vect, zS, fSPhi, fDPhi);
  fGenericPolycone.SetOuterFaceSupport(fR, fZ);
  ComputeBBox();
}

VECCORE_ATT_HOST_DEVICE
void UnplacedGenericPolycone::Extent(Vector3D<Precision> &aMin, Vector3D<Precision> &aMax) const
{
  /*  Algo : Extent along Z direction will be Z[0] and Z[max]
   *         Extent along X and Y will be from the Extent of
   *         cone/tube formed with Maximum outer radius
   */
  Precision rmax = 0.;
  for (unsigned int i = 0; i < fR.size(); i++) {
    if (fR[i] > rmax) {
      rmax = fR[i];
    }
  }
  /* Using the simplest extent for the time being, because
   * currently it's not possible to generate a Cone using
   * factory on GPU. Once it becomes possible then the code
   * written below can be used.
   */
  aMin.Set(-rmax, -rmax, fAMin.z());
  aMax.Set(rmax, rmax, fAMax.z());
  /*auto coneUnplaced = GeoManager::MakeInstance<UnplacedCone>(0., rmax, 0., rmax, 1., fSPhi, fDPhi);
  coneUnplaced->Extent(aMin, aMax);
  aMin.z() = fAMin.z();
  aMax.z() = fAMax.z();*/
}

VECCORE_ATT_HOST_DEVICE
bool UnplacedGenericPolycone::Normal(Vector3D<Precision> const &point, Vector3D<Precision> &norm) const
{
  bool valid = true;
  int index  = fGenericPolycone.GetSectionIndex(point.z() - kTolerance);

  if (index < 0) {
    // Match Polycone: global z end-cap normals are valid surface normals.
    if (index == -1) norm = Vector3D<Precision>(0., 0., -1.);
    if (index == -2) norm = Vector3D<Precision>(0., 0., 1.);
    return valid;
  }
  GenericPolyconeSection const &sec = fGenericPolycone.GetSection(index);
  // At a section seam, only original contour edges are physical boundaries.
  // Adding normalized cone normals cannot cancel their artificial z caps.
  if (size_t(index + 1) < fGenericPolycone.fSections.size() &&
      Abs(point.z() - fGenericPolycone.fZs[index + 1]) < kTolerance) {
    auto const &cone  = *sec.fCoaxialCones->fConeStructVector[0];
    auto const &wedge = cone.fPhiWedge;
    norm.Set(0., 0., 0.);
    if (!cone.IsFullPhi() && !wedge.ContainsWithBoundary(point)) return false;

    const Precision rho = point.Perp();
    Precision area = 0., radial = 0., axial = 0.;
    Precision previousRadial = 0., previousAxial = 0.;
    bool insideRZ = false, onRZ = false;
    for (int i = 0, j = fNumRZ - 1; i < fNumRZ; j = i++) {
      const Precision dr = fR[i] - fR[j], dz = fZ[i] - fZ[j];
      const Precision rr = rho - fR[j], zz = point.z() - fZ[j];
      // Signed area determines contour winding; translation reduces cancellation.
      area += (fR[j] - fR[0]) * (fZ[i] - fZ[0]) - (fR[i] - fR[0]) * (fZ[j] - fZ[0]);
      const Precision cross = rr * dz - zz * dr;
      // Half-open crossings locate phi-only surface points in the RZ polygon.
      if (!cone.IsFullPhi() && ((fZ[j] > point.z()) != (fZ[i] > point.z())) && (dz > 0. ? cross < 0. : cross > 0.))
        insideRZ = !insideRZ;

      const Precision length2 = dr * dr + dz * dz;
      if (length2 == 0.) continue;
      const Precision tolerance = dz == 0. ? kHalfTolerance : kConeTolerance;
      const Precision along     = rr * dr + zz * dz;
      // Squared distance to the finite segment, without a division per edge.
      if (along < 0.) {
        if (rr * rr + zz * zz >= tolerance * tolerance) continue;
      } else if (along > length2) {
        if ((rr - dr) * (rr - dr) + (zz - dz) * (zz - dz) >= tolerance * tolerance) continue;
      } else if (cross * cross >= tolerance * tolerance * length2) {
        continue;
      }
      onRZ = true;
      // A contour segment on the rotation axis is not an exposed lateral face.
      if (fR[i] == 0. && fR[j] == 0.) continue;
      const Precision inverseLength = 1. / Sqrt(length2);
      const Precision nr = dz * inverseLength, nz = -dr * inverseLength;
      // Collinear incident edges describe one face, including at a phi edge.
      if (Abs(nr - previousRadial) + Abs(nz - previousAxial) > 16. * kEpsilon) {
        radial += nr;
        axial += nz;
        previousRadial = nr;
        previousAxial  = nz;
      }
    }
    if (area < 0.) {
      radial = -radial;
      axial  = -axial;
    }
    const Precision inverseRho = 1. / NonZero(rho);
    norm.Set(radial * point.x() * inverseRho, radial * point.y() * inverseRho, axial);
    if (!cone.IsFullPhi() && (insideRZ || onRZ)) {
      if (wedge.IsOnSurface1(point)) norm -= wedge.GetNormal1();
      if (wedge.IsOnSurface2(point)) norm -= wedge.GetNormal2();
    }
    const Precision magnitude2 = norm.Mag2();
    if (magnitude2 == 0.) return false;
    norm /= Sqrt(magnitude2);
    return true;
  }
  // CoaxialCones already normalizes the accumulated valid surface normal.
  return sec.fCoaxialCones->Normal(point - Vector3D<Precision>(0, 0, sec.fShift), norm);
}

Vector3D<Precision> UnplacedGenericPolycone::SamplePointOnSurface() const
{
  /* Algo : First select the section
   *        From the selected section select the cone
   *        Sample the point from selected cone and return the sampled point
   */
  for (int attempt = 0; attempt < 1024; ++attempt) {
    int sectionSelection                        = (int)RNG::Instance().uniform(0., fGenericPolycone.GetNSections());
    const GenericPolyconeSection &section       = fGenericPolycone.GetSection(sectionSelection);
    CoaxialConesStruct<Precision> *coaxialCones = section.fCoaxialCones;
    Vector<ConeStruct<Precision> *> coneStructVector = coaxialCones->fConeStructVector;
    int coneSelection                                = (int)RNG::Instance().uniform(0., coneStructVector.size());
    ConeStruct<Precision> *coneStruct                = coneStructVector[coneSelection];
    SUnplacedCone<ConeTypes::UniversalCone> coneUnplaced(coneStruct->fRmin1, coneStruct->fRmax1, coneStruct->fRmin2,
                                                         coneStruct->fRmax2, coneStruct->fDz, coneStruct->fSPhi,
                                                         coneStruct->fDPhi);
    auto sample = coneUnplaced.SamplePointOnSurface();
    sample.z() += section.fShift;

    // Section-local cone surfaces may be internal shared z planes; the public
    // sampler contract requires points on the full generic-polycone surface.
    if (Inside(sample) == vecgeom::EInside::kSurface) return sample;
  }

  return VUnplacedVolume::SamplePointOnSurface();
}

#ifndef VECCORE_CUDA
bool UnplacedGenericPolycone::CalculateClippedExtent(Axis axis, VoxelLimits const &limits,
                                                     Transformation3D const &transform, Precision &min,
                                                     Precision &max) const
{
  BoundingEnvelope envelope;
  if (!AddGenericPolyconeEnvelope(envelope, *this))
    return VUnplacedVolume::CalculateClippedExtent(axis, limits, transform, min, max);
  return envelope.CalculateExtent(axis, limits, transform, min, max);
}
#endif

VECCORE_ATT_HOST_DEVICE
void UnplacedGenericPolycone::Print() const
{
  // Provided Elliptical Cone Parameters as done for Tube below
  // printf("GenericPolycone {%.2f, %.2f, %.2f}", fGenericPolycone.fDx, fGenericPolycone.fDy, fGenericPolycone.fDz);
  printf("UnplacedGenericPolycone:  {StartPhi: %f   DeltaPhi: %f   N_RZpoints: %d}\n", fSPhi, fDPhi, fNumRZ);
  printf("        -------------------------------------------\n");
  for (int i = 0; i < fNumRZ; i++) {
    printf("        RZpoint #%i:  %f, %f\n", i, fR[i], fZ[i]);
  }
  printf("        -------------------------------------------\n");
}

void UnplacedGenericPolycone::Print(std::ostream &os) const
{
  // Provided Elliptical Cone Parameters as done for Tube below
  // os << "GenericPolycone {" << fEllipticalTube.fDx << ", " << fEllipticalTube.fDy << ", " << fEllipticalTube.fDz <<
  // "}";
  os << "---------------------------------------------------\n"
     << " Solid type: GenericPolycone\n"
     << " Parameters: \n"
     << "    starting phi angle : " << fSPhi << " radians \n"
     << "    ending phi angle   : " << fDPhi << " radians \n";

  os << "    number of RZ points: " << fNumRZ << "\n"
     << "              RZ values (corners): \n";
  for (int i = 0; i < fNumRZ; i++) {
    os << "                         " << fR[i] << ", " << fZ[i] << "\n";
  }
  os << "-----------------------------------------------------------\n";
}

#ifndef VECCORE_CUDA
VPlacedVolume *UnplacedGenericPolycone::Create(LogicalVolume const *const logical_volume,
                                               Transformation3D const *const transformation,
                                               VPlacedVolume *const placement)
{
  if (placement) {
    new (placement) SpecializedGenericPolycone(logical_volume, transformation);
    return placement;
  }
  return new SpecializedGenericPolycone(logical_volume, transformation);
}

VPlacedVolume *UnplacedGenericPolycone::SpecializedVolume(LogicalVolume const *const volume,
                                                          Transformation3D const *const transformation,
                                                          VPlacedVolume *const placement) const
{
  return VolumeFactory::CreateByTransformation<UnplacedGenericPolycone>(volume, transformation, placement);
}
#else

VECCORE_ATT_DEVICE
VPlacedVolume *UnplacedGenericPolycone::Create(LogicalVolume const *const logical_volume,
                                               Transformation3D const *const transformation, const int id,
                                               const int copy_no, const int child_id, VPlacedVolume *const placement)
{
  if (placement) {
    new (placement) SpecializedGenericPolycone(logical_volume, transformation, id, copy_no, child_id);
    return placement;
  }
  return new SpecializedGenericPolycone(logical_volume, transformation, id, copy_no, child_id);
}

VECCORE_ATT_DEVICE
VPlacedVolume *UnplacedGenericPolycone::SpecializedVolume(LogicalVolume const *const volume,
                                                          Transformation3D const *const transformation, const int id,
                                                          const int copy_no, const int child_id,
                                                          VPlacedVolume *const placement) const
{
  return VolumeFactory::CreateByTransformation<UnplacedGenericPolycone>(volume, transformation, id, copy_no, child_id,
                                                                        placement);
}

#endif

#ifdef VECGEOM_CUDA_INTERFACE

DevicePtr<cuda::VUnplacedVolume> UnplacedGenericPolycone::CopyToGpu() const
{
  return CopyToGpuImpl<UnplacedGenericPolycone>();
}

DevicePtr<cuda::VUnplacedVolume> UnplacedGenericPolycone::CopyToGpu(
    DevicePtr<cuda::VUnplacedVolume> const in_gpu_ptr) const
{
  Precision *z_gpu_ptr = AllocateOnGpu<Precision>(fNumRZ * sizeof(Precision));
  Precision *r_gpu_ptr = AllocateOnGpu<Precision>(fNumRZ * sizeof(Precision));

  vecgeom::CopyToGpu(&fZ[0], z_gpu_ptr, sizeof(Precision) * fNumRZ);
  vecgeom::CopyToGpu(&fR[0], r_gpu_ptr, sizeof(Precision) * fNumRZ);

  DevicePtr<cuda::VUnplacedVolume> gpuGenericPolycone =
      CopyToGpuImpl<UnplacedGenericPolycone>(in_gpu_ptr, fSPhi, fDPhi, fNumRZ, r_gpu_ptr, z_gpu_ptr);
  // remove temporary space from GPU
  FreeFromGpu(z_gpu_ptr);
  FreeFromGpu(r_gpu_ptr);

  return gpuGenericPolycone;
}

#endif // VECGEOM_CUDA_INTERFACE

} // namespace VECGEOM_IMPL_NAMESPACE

#ifdef VECCORE_CUDA

namespace cxx {

template size_t DevicePtr<cuda::UnplacedGenericPolycone>::SizeOf();
template void DevicePtr<cuda::UnplacedGenericPolycone>::Construct(Precision phiStart, Precision phiTotal, int numRZ,
                                                                  Precision *r, Precision *z) const;

} // namespace cxx

#endif
} // namespace vecgeom
