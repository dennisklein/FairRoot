/// \file UnplacedRootVolume.cpp
/// \author Johannes de Fine Licht (johannes.definelicht@cern.ch)

#include "UnplacedRootVolume.h"
#include "PlacedRootVolume.h"
#include <stdio.h>

namespace vecgeom::tgeo2vecgeom {

void UnplacedRootVolume::Print() const
{
  printf("UnplacedRootVolume");
}

void UnplacedRootVolume::Print(std::ostream &os) const
{
  os << "UnplacedRootVolume";
}

VPlacedVolume *UnplacedRootVolume::SpecializedVolume(LogicalVolume const *const volume,
                                                     Transformation3D const *const transformation,
#if !TGEO2VECGEOM_VECGEOM_V2
                                                     const TranslationCode trans_code, const RotationCode rot_code,
#endif
                                                     VPlacedVolume *const placement) const
{
  if (placement) {
    return new (placement) PlacedRootVolume(fRootShape, volume, transformation);
  }
  return new PlacedRootVolume(fRootShape, volume, transformation);
}

} // namespace vecgeom::tgeo2vecgeom
