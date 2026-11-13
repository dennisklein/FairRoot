/// \file VecGeomCompat.h
///
/// Single place where the VecGeom 1.x / 2.x API split is named.
///
/// VecGeom 2.0 reworked the volume interfaces this package derives from:
///  * the SOA3D and SIMD (`*Vec`) entry points of VPlacedVolume are gone, and
///    DistanceToIn/DistanceToOut/PlacedDistanceToOut gained an overload taking a
///    SurfaceHitView out-parameter;
///  * VUnplacedVolume::SpecializedVolume() lost its TranslationCode/RotationCode
///    arguments (there is no compile-time specialisation any more), and Normal()
///    and Extent() became pure virtual;
///  * LogicalVolume lost the opaque fMaterialPtr slot;
///  * Transformation3D::Inverse() returns the inverse instead of filling an
///    out-parameter, and NavigationState took over Index2PVolumeConverter.
///
/// Rather than sprinkling raw version tests through the sources, every guard
/// below is written as `#if TGEO2VECGEOM_VECGEOM_V2`, so that dropping VecGeom
/// 1.x support later is a matter of deleting the `#else` branches.

#pragma once

#include "VecGeom/base/Version.h"

/// VECGEOM_VERSION is (major << 16) | (minor << 8) | patch.
#if VECGEOM_VERSION >= 0x020000
#define TGEO2VECGEOM_VECGEOM_V2 1
#else
#define TGEO2VECGEOM_VECGEOM_V2 0
#endif
