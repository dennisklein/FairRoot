/// @file UnplacedExtruded.cpp
/// @author Mihaela Gheata (mihaela.gheata@cern.ch)

#include "VecGeom/volumes/Tessellated.h"
#include "VecGeom/base/Triangulation.h"
#include "VecGeom/volumes/UnplacedExtruded.h"
#include "VecGeom/volumes/UnplacedSExtruVolume.h"
#include "VecGeom/volumes/PlanarPolygonUtilities.h"
#include "VecGeom/volumes/SpecializedExtruded.h"
#include "VecGeom/volumes/utilities/GenerationUtilities.h"
#include "VecGeom/base/RNG.h"

#include "VecGeom/management/VolumeFactory.h"
#include <array>
#include <cstddef>
#include <vector>

#ifndef VECCORE_CUDA
#include "VecGeom/volumes/UnplacedImplAs.h"
#include "VecGeom/volumes/utilities/BoundingEnvelope.h"
#endif

#ifndef VECCORE_CUDA
#ifdef VECGEOM_ROOT
#include <TGeoXtru.h>
#endif

#ifdef VECGEOM_GEANT4
#include <G4ExtrudedSolid.hh>
#include <G4TessellatedSolid.hh>
#include <G4TriangularFacet.hh>
#endif
#endif

namespace vecgeom {

inline namespace VECGEOM_IMPL_NAMESPACE {

#ifndef VECCORE_CUDA
namespace {

/// @brief Triangulate a PlanarPolygon with the extruded ear-clipping loop.
template <typename TriangleList>
bool TriangulatePlanarPolygon(PlanarPolygon const &polygon, TriangleList &triangles)
{
  return TriangulateEarClipping(
      polygon.GetNVertices(), triangles,
      [](size_t i1, size_t i2, size_t i3) { return std::array<size_t, 3>{{i1, i2, i3}}; },
      [&polygon](size_t first, size_t second, size_t third) {
        return PlanarPolygonUtilities::IsConvexCorner(polygon, first, second, third);
      },
      [&polygon](size_t first, size_t second, size_t third, size_t index) {
        return PlanarPolygonUtilities::IsPointInsideTriangle(polygon, first, second, third, index);
      });
}

/// @brief Map a base polygon vertex to one extrusion section.
Vector3D<Precision> VertexAtSection(PlanarPolygon const &polygon, size_t vertex, XtruSection const &section)
{
  return section.fOrigin + section.fScale * PlanarPolygonUtilities::Vertex2D(polygon, vertex);
}

/// @brief Add triangular-prism envelope pieces for all adjacent section pairs.
bool AddExtrudedPolygonEnvelope(BoundingEnvelope &envelope, PlanarPolygon const &polygon,
                                std::vector<XtruSection> const &sections)
{
  if (polygon.GetNVertices() < 3 || sections.size() < 2) return false;

  for (auto const &section : sections) {
    if (section.fScale <= 0.) return false;
  }

  std::vector<std::array<size_t, 3>> triangles;
  if (!TriangulatePlanarPolygon(polygon, triangles)) return false;

  // Match Geant4's base-triangle decomposition, but add one prism per adjacent
  // section pair so each lateral face follows the piecewise-linear extrusion.
  bool added_piece = false;
  for (auto const &triangle : triangles) {
    for (size_t isect = 0; isect + 1 < sections.size(); ++isect) {
      if (sections[isect + 1].fOrigin.z() < sections[isect].fOrigin.z()) return false;
      if (sections[isect + 1].fOrigin.z() < sections[isect].fOrigin.z() + kTolerance) continue;

      Vector3D<Precision> lower[3];
      Vector3D<Precision> upper[3];
      for (int vertex = 0; vertex < 3; ++vertex) {
        lower[vertex] = VertexAtSection(polygon, triangle[vertex], sections[isect]);
        upper[vertex] = VertexAtSection(polygon, triangle[vertex], sections[isect + 1]);
      }
      envelope.AddLocalPrism(lower, 3, upper, 3);
      added_piece = true;
    }
  }
  return added_piece;
}

/// @brief Add envelope pieces for ExtrudedStruct using its existing mesh triangulation.
bool AddExtrudedStructEnvelope(BoundingEnvelope &envelope, ExtrudedStruct const &xtru)
{
  const size_t nvertices = xtru.GetNVertices();
  const size_t nsections = xtru.GetNSections();
  if (nvertices < 3 || nsections < 2) return false;

  for (size_t isect = 0; isect < nsections; ++isect) {
    if (xtru.GetSection(static_cast<int>(isect)).fScale <= 0.) return false;
  }

  const auto mesh_helper  = xtru.GetMeshHelper();
  const size_t ntriangles = nvertices - 2;
  bool added_piece        = false;
  for (size_t itri = 0; itri < ntriangles; ++itri) {
    size_t mesh_indices[3] = {};
    mesh_helper.GetFacetVertices(itri, mesh_indices);
    size_t triangle[3] = {mesh_indices[0] % nvertices, mesh_indices[1] % nvertices, mesh_indices[2] % nvertices};

    for (size_t isect = 0; isect + 1 < nsections; ++isect) {
      const auto lower_section = xtru.GetSection(static_cast<int>(isect));
      const auto upper_section = xtru.GetSection(static_cast<int>(isect + 1));
      if (upper_section.fOrigin.z() < lower_section.fOrigin.z()) return false;
      if (upper_section.fOrigin.z() < lower_section.fOrigin.z() + kTolerance) continue;

      Vector3D<Precision> lower[3];
      Vector3D<Precision> upper[3];
      for (int vertex = 0; vertex < 3; ++vertex) {
        lower[vertex] = xtru.VertexToSection(triangle[vertex], isect);
        upper[vertex] = xtru.VertexToSection(triangle[vertex], isect + 1);
      }
      envelope.AddLocalPrism(lower, 3, upper, 3);
      added_piece = true;
    }
  }
  return added_piece;
}

/// @brief Adapt a two-plane PolygonalShell to the section format used here.
std::vector<XtruSection> MakeSExtruSections(PolygonalShell const &shell)
{
  std::vector<XtruSection> sections(2);
  sections[0].fOrigin.Set(0., 0., shell.GetLowerZ());
  sections[0].fScale = 1.;
  sections[1].fOrigin.Set(0., 0., shell.GetUpperZ());
  sections[1].fScale = 1.;
  return sections;
}

/// @brief Add envelope pieces for a simple SExtru polygonal shell.
bool AddPolygonShellEnvelope(BoundingEnvelope &envelope, PolygonalShell const &shell)
{
  return AddExtrudedPolygonEnvelope(envelope, shell.GetPolygon(), MakeSExtruSections(shell));
}

} // namespace

#ifdef VECGEOM_ROOT
TGeoShape const *UnplacedExtruded::ConvertToRoot(char const *label) const
{
  size_t nvert = GetNVertices();
  size_t nsect = GetNSections();

  // if(nsect > 1){
  double *x = new double[nvert];
  double *y = new double[nvert];
  for (size_t i = 0; i < nvert; ++i) {
    Precision xcrt, ycrt;
    GetVertex(i, xcrt, ycrt);
    x[i] = xcrt;
    y[i] = ycrt;
  }
  TGeoXtru *xtru = new TGeoXtru(nsect);
  xtru->DefinePolygon(nvert, x, y);
  for (size_t i = 0; i < nsect; ++i) {
    XtruSection sect = GetSection(i);
    xtru->DefineSection(i, sect.fOrigin.z(), sect.fOrigin.x(), sect.fOrigin.y(), sect.fScale);
  }
  delete[] x;
  delete[] y;
  return xtru;
}
#endif

#ifdef VECGEOM_GEANT4
G4VSolid const *UnplacedExtruded::ConvertToGeant4(char const *label) const
{
  std::vector<G4TwoVector> polygon;
  Precision x, y;
  size_t nvert = GetNVertices();
  for (size_t i = 0; i < nvert; ++i) {
    GetVertex(i, x, y);
    polygon.push_back(G4TwoVector(x, y));
  }
  std::vector<G4ExtrudedSolid::ZSection> sections;
  size_t nsect = GetNSections();
  for (size_t i = 0; i < nsect; ++i) {
    XtruSection sect = GetSection(i);
    sections.push_back(
        G4ExtrudedSolid::ZSection(sect.fOrigin.z(), G4TwoVector(sect.fOrigin.x(), sect.fOrigin.y()), sect.fScale));
  }
  G4ExtrudedSolid *g4xtru = new G4ExtrudedSolid(label, polygon, sections);
  return g4xtru;
}
#endif
#endif

template <>
UnplacedExtruded *Maker<UnplacedExtruded>::MakeInstance(const size_t nvertices, XtruVertex2 const *vertices,
                                                        const int nsections, XtruSection const *sections)
{

#ifndef VECGEOM_NO_SPECIALIZATION
  bool isSExtru = false;
  for (int i = 0; i < (nsections - 1); i++) {
    if (i == 0) {
      isSExtru = ((sections[i].fOrigin - sections[i + 1].fOrigin).Perp2() < kTolerance &&
                  Abs(sections[i].fScale - sections[i + 1].fScale) < kTolerance);
    } else {
      isSExtru &= ((sections[i].fOrigin - sections[i + 1].fOrigin).Perp2() < kTolerance &&
                   Abs(sections[i].fScale - sections[i + 1].fScale) < kTolerance);
    }
    if (!isSExtru) break;
  }
  if (isSExtru) {
    Precision *x = new Precision[nvertices];
    Precision *y = new Precision[nvertices];
    for (size_t i = 0; i < nvertices; ++i) {
      x[i] = vertices[i].x;
      y[i] = vertices[i].y;
    }
    Precision zmin = sections[0].fOrigin.z();
    Precision zmax = sections[nsections - 1].fOrigin.z();
    return new SUnplacedImplAs<UnplacedExtruded, UnplacedSExtruVolume>(nvertices, x, y, zmin, zmax);
  } else {
    return new UnplacedExtruded(nvertices, vertices, nsections, sections);
  }
#else
  return new UnplacedExtruded(nvertices, vertices, nsections, sections);
#endif
}

void UnplacedExtruded::Print() const
{
#ifdef VECCORE_CUDA
  printf("UnplacedExtruded");
#else
  std::cerr << "UnplacedExtruded: vertices {";
  int nvert = GetNVertices();
  Precision x, y;
  for (int i = 0; i < nvert - 1; ++i) {
    GetVertex(i, x, y);
    std::cerr << "(" << x << ", " << y << "), ";
  }
  GetVertex(nvert - 1, x, y);
  std::cerr << "(" << x << ", " << y << ")}\n";
  std::cerr << "sections:\n";
  int nsect = GetNSections();
  for (int i = 0; i < nsect; ++i) {
    XtruSection sect = GetSection(i);
    std::cerr << "orig: (" << sect.fOrigin.x() << ", " << sect.fOrigin.y() << ", " << sect.fOrigin.z()
              << ") scl = " << sect.fScale << std::endl;
  }
#endif
}

void UnplacedExtruded::Print(std::ostream &os) const
{
  os << "UnplacedExtruded: vertices {";
  int nvert = GetNVertices();
  Precision x, y;
  for (int i = 0; i < nvert - 1; ++i) {
    GetVertex(i, x, y);
    os << "(" << x << ", " << y << "), ";
  }
  GetVertex(nvert - 1, x, y);
  os << "(" << x << ", " << y << ")}\n";
  os << "sections:\n";
  int nsect = GetNSections();
  for (int i = 0; i < nsect; ++i) {
    XtruSection sect = GetSection(i);
    os << "orig: (" << sect.fOrigin.x() << ", " << sect.fOrigin.y() << ", " << sect.fOrigin.z()
       << ") scl = " << sect.fScale << std::endl;
  }
}

void UnplacedExtruded::Extent(Vector3D<Precision> &aMin, Vector3D<Precision> &aMax) const
{
  // TODO: this differentiated treatment seems not to be needed
  if (fXtru.fIsSxtru) {
    fXtru.fSxtruHelper.Extent(aMin, aMax);
  } else {
    aMin = fXtru.fTslRuntimeHelper.fMinExtent;
    aMax = fXtru.fTslRuntimeHelper.fMaxExtent;
  }
}

#ifndef VECCORE_CUDA
bool UnplacedExtruded::CalculateClippedExtent(Axis axis, VoxelLimits const &limits, Transformation3D const &transform,
                                              Precision &min, Precision &max) const
{
  if (fXtru.fIsSxtru) {
    BoundingEnvelope envelope;
    if (!AddPolygonShellEnvelope(envelope, fXtru.fSxtruHelper))
      return VUnplacedVolume::CalculateClippedExtent(axis, limits, transform, min, max);
    return envelope.CalculateExtent(axis, limits, transform, min, max);
  }

  BoundingEnvelope envelope;
  if (!AddExtrudedStructEnvelope(envelope, fXtru))
    return VUnplacedVolume::CalculateClippedExtent(axis, limits, transform, min, max);

  return envelope.CalculateExtent(axis, limits, transform, min, max);
}

bool UnplacedSExtruVolume::CalculateClippedExtent(Axis axis, VoxelLimits const &limits,
                                                  Transformation3D const &transform, Precision &min,
                                                  Precision &max) const
{
  BoundingEnvelope envelope;
  if (!AddPolygonShellEnvelope(envelope, fPolyShell))
    return VUnplacedVolume::CalculateClippedExtent(axis, limits, transform, min, max);
  return envelope.CalculateExtent(axis, limits, transform, min, max);
}
#endif

Precision UnplacedExtruded::Capacity() const
{
  if (fXtru.fCubicVolume != 0.) return fXtru.fCubicVolume;

  if (fXtru.fIsSxtru) {
    fXtru.fCubicVolume =
        fXtru.fSxtruHelper.GetPolygon().Area() * (fXtru.fSxtruHelper.GetUpperZ() - fXtru.fSxtruHelper.GetLowerZ());
  } else {
    int size = fXtru.fTslRuntimeHelper.fNFacets;
    for (int i = 0; i < size; ++i) {
      const auto &facet = fXtru.fTslRuntimeHelper.fFacets[i];
      Precision area    = facet.SurfaceArea();
      fXtru.fCubicVolume += area * (facet.fVertices[0].Dot(facet.fNormal));
    }
    fXtru.fCubicVolume /= 3.;
  }
  return fXtru.fCubicVolume;
}

Precision UnplacedExtruded::SurfaceArea() const
{
  if (fXtru.fSurfaceArea != 0.) return fXtru.fSurfaceArea;

  if (fXtru.fIsSxtru) {
    fXtru.fSurfaceArea = fXtru.fSxtruHelper.SurfaceArea() + 2. * fXtru.fSxtruHelper.GetPolygon().Area();
  } else {
    int size = fXtru.fTslRuntimeHelper.fNFacets;
    for (int i = 0; i < size; ++i) {
      const auto &facet = fXtru.fTslRuntimeHelper.fFacets[i];
      fXtru.fSurfaceArea += facet.SurfaceArea();
    }
  }
  return fXtru.fSurfaceArea;
}

int UnplacedExtruded::ChooseSurface() const
{
  int choice       = 0; // 0 = zm, 1 = zp, 2 = ym, 3 = yp, 4 = xm, 5 = xp
  Precision Stotal = SurfaceArea();

  // random value to choose surface to place the point
  Precision rand = RNG::Instance().uniform() * Stotal;

  while (rand > fXtru.fTslRuntimeHelper.fFacets[choice].SurfaceArea())
    rand -= fXtru.fTslRuntimeHelper.fFacets[choice].SurfaceArea(), choice++;

  return choice;
}

Vector3D<Precision> UnplacedExtruded::SamplePointOnSurface() const
{
  int surface  = ChooseSurface();
  Precision r1 = RNG::Instance().uniform(0.0, 1.0);
  Precision r2 = RNG::Instance().uniform(0.0, 1.0);
  if (r1 + r2 > 1.) {
    r1 = 1. - r1;
    r2 = 1. - r2;
  }
  const auto &facet = fXtru.fTslRuntimeHelper.fFacets[surface];
  return (facet.fVertices[0] + r1 * (facet.fVertices[1] - facet.fVertices[0]) +
          r2 * (facet.fVertices[2] - facet.fVertices[0]));
}

bool UnplacedExtruded::Normal(Vector3D<Precision> const &point, Vector3D<Precision> &norm) const
{
  // Redirect to normal implementation
  bool valid = false;
  if (fXtru.fIsSxtru) {
    norm = SExtruImplementation::NormalKernel(fXtru.fSxtruHelper, point, valid);
  } else {
    norm = TessellatedImplementation::NormalKernel<Precision>(fXtru.fTslRuntimeHelper, point, valid);
  }
  return valid;
}

#ifndef VECCORE_CUDA
SolidMesh *UnplacedExtruded::CreateMesh3D(Transformation3D const &trans, size_t nSegments) const
{

  typedef Vector3D<Precision> Vec_t;

  SolidMesh *sm = new SolidMesh();

  size_t n         = GetNVertices();
  size_t nSections = GetNSections();

  Vec_t *vertices = new Vec_t[nSections * (n + 1)];
  size_t idx      = 0;
  for (size_t i = 0; i < nSections; i++) {
    for (size_t j = n; j > 0; j--) {
      vertices[idx++] = GetStruct().VertexToSection(j - 1, i);
    }
    vertices[idx++] = GetStruct().VertexToSection(n - 1, i);
  }

  sm->ResetMesh(nSections * (n + 1), (nSections - 1) * (n) + 2);
  sm->SetVertices(vertices, nSections * (n + 1));
  delete[] vertices;
  sm->TransformVertices(trans);

  std::vector<size_t> indices;
  for (size_t i = n; i > 0; i--) {
    indices.push_back(i - 1);
  }
  sm->AddPolygon(n, indices, GetStruct().IsConvexPolygon());

  indices.clear();

  for (size_t i = 0, k = (nSections - 1) * (n + 1); i < n; i++, k++) {
    indices.push_back(k);
  }

  sm->AddPolygon(n, indices, GetStruct().IsConvexPolygon());

  size_t k = 0;
  for (size_t i = 0; i < nSections - 1; i++, k++) {
    for (size_t j = 0; j < n; j++, k++) {
      sm->AddPolygon(4, {k, k + 1, k + 1 + n + 1, k + n + 1}, true);
    }
  }

  return sm;
}
#endif

#ifdef VECCORE_CUDA
VECCORE_ATT_DEVICE VPlacedVolume *UnplacedExtruded::Create(LogicalVolume const *const logical_volume,
                                                           Transformation3D const *const transformation, const int id,
                                                           const int copy_no, const int child_id,
                                                           VPlacedVolume *const placement)
{
  if (placement) {
    new (placement) SpecializedExtruded(logical_volume, transformation, id, copy_no, child_id);
    return placement;
  }
  return new SpecializedExtruded(logical_volume, transformation, id, copy_no, child_id);
}
#else
VPlacedVolume *UnplacedExtruded::Create(LogicalVolume const *const logical_volume,
                                        Transformation3D const *const transformation, VPlacedVolume *const placement)
{
  if (placement) {
    new (placement) SpecializedExtruded(logical_volume, transformation);
    return placement;
  }
  return new SpecializedExtruded(logical_volume, transformation);
}
#endif

VECCORE_ATT_DEVICE
VPlacedVolume *UnplacedExtruded::SpecializedVolume(LogicalVolume const *const volume,
                                                   Transformation3D const *const transformation,
#ifdef VECCORE_CUDA
                                                   const int id, const int copy_no, const int child_id,
#endif
                                                   VPlacedVolume *const placement) const
{

  return VolumeFactory::CreateByTransformation<UnplacedExtruded>(volume, transformation,
#ifdef VECCORE_CUDA
                                                                 id, copy_no, child_id,
#endif
                                                                 placement);
}

std::ostream &UnplacedExtruded::StreamInfo(std::ostream &os) const
{
  int oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "     *** Dump for solid - " << GetEntityType() << " ***\n"
     << "     ===================================================\n"
     << " Solid type: Extruded\n"
     << " Parameters: \n"
     << "-----------------------------------------------------------\n";
  os.precision(oldprc);
  return os;
}

#ifdef VECGEOM_CUDA_INTERFACE

DevicePtr<cuda::VUnplacedVolume> UnplacedExtruded::CopyToGpu(DevicePtr<cuda::VUnplacedVolume> const in_gpu_ptr) const
{
  if (GetStruct().fIsSxtru) {
    auto const &shell      = GetStruct().fSxtruHelper;
    auto const &vertices   = shell.GetPolygon().GetVertices();
    Precision const *x_cpu = vertices.x();
    Precision const *y_cpu = vertices.y();
    const auto nvertices   = vertices.size();
    Precision *x_gpu_ptr   = AllocateOnGpu<Precision>(nvertices * sizeof(Precision));
    Precision *y_gpu_ptr   = AllocateOnGpu<Precision>(nvertices * sizeof(Precision));
    vecgeom::CopyToGpu(x_cpu, x_gpu_ptr, sizeof(Precision) * nvertices);
    vecgeom::CopyToGpu(y_cpu, y_gpu_ptr, sizeof(Precision) * nvertices);

    auto final_gpu_ptr = CopyToGpuImpl<UnplacedExtruded>(in_gpu_ptr, static_cast<int>(nvertices), x_gpu_ptr, y_gpu_ptr,
                                                         shell.GetLowerZ(), shell.GetUpperZ(), fGlobalConvexity);

    FreeFromGpu(x_gpu_ptr);
    FreeFromGpu(y_gpu_ptr);
    return final_gpu_ptr;
  }

  auto const &runtime = GetStruct().fTslRuntimeHelper;
  VECGEOM_VALIDATE(runtime.fNFacets > 0 && runtime.fFacets != nullptr && runtime.fBVH != nullptr,
                   << "Attempted to copy UnplacedExtruded without initialized tessellated runtime data");

  // The facet hierarchy belongs to this solid rather than the manager's volume-navigation table,
  // so give it an independent device object. Use the host size because cuda::IndexedBVH is only
  // forward-declared here and both namespace variants have the same layout.
  auto gpu_bvh_ptr = AllocateOnGpu<cuda::IndexedBVH<float>>(sizeof(IndexedBVH<float>));
  runtime.fBVH->CopyToGpu(gpu_bvh_ptr);

  size_t nfacets      = runtime.fNFacets;
  auto gpu_facets_ptr = AllocateOnGpu<cuda::TriangularTile<Precision>>(sizeof(TriangularTile<Precision>) * nfacets);
  vecgeom::CopyToGpu((char *)runtime.fFacets, (char *)gpu_facets_ptr, sizeof(TriangularTile<Precision>) * nfacets);

  auto gpu_supports_ptr = AllocateOnGpu<unsigned char>(sizeof(unsigned char) * nfacets);
  vecgeom::CopyToGpu(runtime.fFacetSupports, gpu_supports_ptr, sizeof(unsigned char) * nfacets);

  const bool use_tsl_sections         = GetStruct().fUseTslSections;
  const bool use_projection           = GetStruct().fProjectionSafetySqr > 0.;
  const size_t nsections              = (use_tsl_sections || use_projection) ? GetStruct().fSections.size() : 0;
  const size_t ntsl_sections          = use_tsl_sections ? GetStruct().fTslSections.size() : 0;
  Precision *gpu_zplanes_ptr          = nullptr;
  cuda::XtruSection *gpu_sections_ptr = nullptr;
  cuda::TessellatedSection<Precision> **gpu_tsl_section_ptrs = nullptr;

  if (nsections > 0) {
    gpu_zplanes_ptr = AllocateOnGpu<Precision>(sizeof(Precision) * nsections);
    vecgeom::CopyToGpu(GetStruct().fZPlanes, gpu_zplanes_ptr, sizeof(Precision) * nsections);

    gpu_sections_ptr = AllocateOnGpu<cuda::XtruSection>(sizeof(XtruSection) * nsections);
    vecgeom::CopyToGpu((char *)GetStruct().fSections.begin(), (char *)gpu_sections_ptr,
                       sizeof(XtruSection) * nsections);
  }

  if (use_tsl_sections) {
    VECGEOM_VALIDATE(nsections > 1 && ntsl_sections == nsections - 1,
                     << "Extruded tessellated section helper has inconsistent section counts");

    size_t total_section_facets = 0;
    for (size_t isect = 0; isect < ntsl_sections; ++isect)
      total_section_facets += GetStruct().fTslSections[isect]->GetNfacets();

    std::vector<QuadrilateralFacet<Precision>> section_facets;
    section_facets.reserve(total_section_facets);
    for (size_t isect = 0; isect < ntsl_sections; ++isect) {
      auto const *section = GetStruct().fTslSections[isect];
      auto const *facets  = section->GetFacetData();
      section_facets.insert(section_facets.end(), facets, facets + section->GetNfacets());
    }

    auto gpu_section_facets_ptr =
        AllocateOnGpu<cuda::Tile<4, Precision>>(total_section_facets * sizeof(QuadrilateralFacet<Precision>));
    vecgeom::CopyToGpu((char *)section_facets.data(), (char *)gpu_section_facets_ptr,
                       total_section_facets * sizeof(QuadrilateralFacet<Precision>));

    const size_t section_size = DevicePtr<cuda::TessellatedSection<Precision>>::SizeOf();
    auto gpu_tsl_sections_ptr = AllocateOnGpu<cuda::TessellatedSection<Precision>>(ntsl_sections * section_size);
    std::vector<cuda::TessellatedSection<Precision> *> section_ptrs(ntsl_sections);

    size_t facet_offset = 0;
    for (size_t isect = 0; isect < ntsl_sections; ++isect) {
      auto const *section = GetStruct().fTslSections[isect];
      const size_t nquad  = section->GetNfacets();

      auto *facet_slice = reinterpret_cast<cuda::Tile<4, Precision> *>(
          reinterpret_cast<char *>(gpu_section_facets_ptr) + facet_offset * sizeof(QuadrilateralFacet<Precision>));
      auto *section_ptr = reinterpret_cast<cuda::TessellatedSection<Precision> *>(
          reinterpret_cast<char *>(gpu_tsl_sections_ptr) + isect * section_size);
      section_ptrs[isect] = section_ptr;

      DevicePtr<cuda::TessellatedSection<Precision>> gpu_section(section_ptr);
      // Device sections wrap pre-copied facet slices; no section helper is rebuilt on device.
      gpu_section.Construct(facet_slice, nquad, GetStruct().fZPlanes[isect], GetStruct().fZPlanes[isect + 1]);
      facet_offset += nquad;
    }

    gpu_tsl_section_ptrs = AllocateOnGpu<cuda::TessellatedSection<Precision> *>(
        sizeof(cuda::TessellatedSection<Precision> *) * ntsl_sections);
    vecgeom::CopyToGpu((char *)section_ptrs.data(), (char *)gpu_tsl_section_ptrs,
                       sizeof(cuda::TessellatedSection<Precision> *) * ntsl_sections);
  }

  auto final_gpu_ptr = CopyToGpuImpl<UnplacedExtruded>(
      in_gpu_ptr, nfacets, gpu_facets_ptr, gpu_bvh_ptr, fGlobalConvexity, gpu_zplanes_ptr, gpu_sections_ptr, nsections,
      gpu_tsl_section_ptrs, ntsl_sections, use_tsl_sections, gpu_supports_ptr);

  if (use_projection) {
    const auto &vertices   = GetStruct().fPolygon.GetVertices();
    const int polygon_size = GetStruct().fPolygon.GetNVertices();
    auto *polygon_x        = AllocateOnGpu<Precision>(sizeof(Precision) * polygon_size);
    auto *polygon_y        = AllocateOnGpu<Precision>(sizeof(Precision) * polygon_size);
    vecgeom::CopyToGpu(vertices.x(), polygon_x, sizeof(Precision) * polygon_size);
    vecgeom::CopyToGpu(vertices.y(), polygon_y, sizeof(Precision) * polygon_size);
    // Reuse the polygon's existing device constructor in its empty member slot.
    // The constructor owns its copies, so temporary coordinate buffers can be freed.
    const auto *host_bytes    = reinterpret_cast<const char *>(this);
    auto *device_bytes        = reinterpret_cast<char *>(in_gpu_ptr.GetPtr());
    const auto polygon_offset = reinterpret_cast<const char *>(&GetStruct().fPolygon) - host_bytes;
    DevicePtr<cuda::PlanarPolygon> gpu_polygon(static_cast<void *>(device_bytes + polygon_offset));
    gpu_polygon.Construct(polygon_size, polygon_x, polygon_y);
    const auto safety_offset = reinterpret_cast<const char *>(&GetStruct().fProjectionSafetySqr) - host_bytes;
    const auto z_offset      = reinterpret_cast<const char *>(&GetStruct().fProjectionZTolerance) - host_bytes;
    vecgeom::CopyToGpu(host_bytes + safety_offset, device_bytes + safety_offset, sizeof(Precision));
    vecgeom::CopyToGpu(host_bytes + z_offset, device_bytes + z_offset, sizeof(Precision));
    FreeFromGpu(polygon_x);
    FreeFromGpu(polygon_y);
  }

  const auto *this_bytes                  = reinterpret_cast<const char *>(this);
  const auto *runtime_bytes               = reinterpret_cast<const char *>(&runtime);
  const std::size_t runtime_struct_offset = static_cast<std::size_t>(runtime_bytes - this_bytes);
  // Preserve the device facet/BVH pointers set by the GPU constructor; copy only host-computed metadata.
  constexpr std::size_t first_pointer_offset_in_runtime = offsetof(TessellatedRuntimeStruct<Precision>, fFacets);
  vecgeom::CopyToGpu((char *)this + runtime_struct_offset, (char *)(in_gpu_ptr.GetPtr()) + runtime_struct_offset,
                     first_pointer_offset_in_runtime);

  return final_gpu_ptr;
}

DevicePtr<cuda::VUnplacedVolume> UnplacedExtruded::CopyToGpu() const { return CopyToGpuImpl<UnplacedExtruded>(); }

#endif // VECGEOM_CUDA_INTERFACE

} // namespace VECGEOM_IMPL_NAMESPACE

#ifdef VECCORE_CUDA

namespace cxx {

template size_t DevicePtr<cuda::UnplacedExtruded>::SizeOf();
template void DevicePtr<cuda::UnplacedExtruded>::Construct() const;
template void DevicePtr<cuda::UnplacedExtruded>::Construct(size_t, TriangularTile<Precision> *, IndexedBVH<float> *,
                                                           bool) const;
template void DevicePtr<cuda::UnplacedExtruded>::Construct(size_t, TriangularTile<Precision> *, IndexedBVH<float> *,
                                                           bool, Precision *, cuda::XtruSection *, size_t,
                                                           cuda::TessellatedSection<Precision> **, size_t, bool,
                                                           unsigned char *) const;
template void DevicePtr<cuda::PlanarPolygon>::Construct(int, Precision *, Precision *) const;
template void DevicePtr<cuda::UnplacedExtruded>::Construct(int, Precision *, Precision *, Precision, Precision,
                                                           bool) const;

} // namespace cxx

#endif

} // namespace vecgeom
