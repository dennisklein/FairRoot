/// \file TessellatedHelpers.cpp
/// \author Mihaela Gheata (mihaela.gheata@cern.ch)

#include <ostream>
#include "VecGeom/volumes/TessellatedSection.h"

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

#ifndef VECCORE_CUDA
std::ostream &operator<<(std::ostream &os, vecgeom::TessellatedSection<double> const &ts)
{
  os << " tessellated section with " << ts.GetNfacets() << " facets:\n";
  for (size_t i = 0; i < ts.GetNfacets(); ++i) {
    os << ts.GetFacet(i) << "\n";
  }
  return os;
}
#endif
} // namespace VECGEOM_IMPL_NAMESPACE

#ifdef VECCORE_CUDA
namespace cxx {
template size_t DevicePtr<cuda::TessellatedSection<Precision>>::SizeOf();
template void DevicePtr<cuda::TessellatedSection<Precision>>::Construct(cuda::Tile<4, Precision> *, size_t, Precision,
                                                                        Precision) const;
} // namespace cxx
#endif
} // namespace vecgeom
