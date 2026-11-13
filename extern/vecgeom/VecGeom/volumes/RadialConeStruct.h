/// @file RadialConeStruct.h
/// @brief Radial and z data for complete cone sections without phi state.

#ifndef VECGEOM_VOLUMES_RADIALCONESTRUCT_H_
#define VECGEOM_VOLUMES_RADIALCONESTRUCT_H_

#include "VecGeom/base/Global.h"

namespace vecgeom {

VECGEOM_DEVICE_DECLARE_CONV_TEMPLATE(struct, RadialConeStruct, typename);

inline namespace VECGEOM_IMPL_NAMESPACE {

template <typename T = double>
struct RadialConeStruct {
  T fRmin1{0.};
  T fRmax1{0.};
  T fRmin2{0.};
  T fRmax2{0.};
  T fDz{0.};

  Precision fInnerSlope{0.};
  Precision fOuterSlope{0.};
  Precision fInnerOffset{0.};
  Precision fOuterOffset{0.};
  Precision fInnerTolerance{0.};
  Precision fOuterTolerance{0.};
  Precision fSqRmin1{0.}, fSqRmin2{0.};
  Precision fSqRmax1{0.}, fSqRmax2{0.};
  Precision fTolIz{0.}, fTolOz{0.};
  Precision fInnerConeApex{0.};
  Precision fTanInnerApexAngle{0.};
  Precision fOuterConeApex{0.};
  Precision fTanOuterApexAngle{0.};
  Precision fSecRMin{0.};
  Precision fSecRMax{0.};
  Precision fInvSecRMin{0.};
  Precision fInvSecRMax{0.};
  Precision fTanRMin{0.};
  Precision fTanRMax{0.};
  Precision fZNormInner{0.};
  Precision fZNormOuter{0.};

  VECCORE_ATT_HOST_DEVICE
  Precision Capacity(Precision dphi) const
  {
    return (fDz * dphi / 3.) *
           (fRmax1 * fRmax1 + fRmax2 * fRmax2 + fRmax1 * fRmax2 - fRmin1 * fRmin1 - fRmin2 * fRmin2 - fRmin1 * fRmin2);
  }

  VECCORE_ATT_HOST_DEVICE
  void CalculateCached()
  {

    fSqRmin1 = fRmin1 * fRmin1;
    fSqRmax1 = fRmax1 * fRmax1;
    fSqRmin2 = fRmin2 * fRmin2;
    fSqRmax2 = fRmax2 * fRmax2;

    fTanRMin    = (fRmin2 - fRmin1) * 0.5 / fDz;
    fSecRMin    = Sqrt(1.0 + fTanRMin * fTanRMin);
    fInvSecRMin = 1. / NonZero(fSecRMin);
    fTanRMax    = (fRmax2 - fRmax1) * 0.5 / fDz;
    fSecRMax    = Sqrt(1.0 + fTanRMax * fTanRMax);
    fInvSecRMax = 1. / NonZero(fSecRMax);

    fInnerSlope     = -(fRmin1 - fRmin2) / (2. * fDz);
    fOuterSlope     = -(fRmax1 - fRmax2) / (2. * fDz);
    fInnerOffset    = fRmin2 - fInnerSlope * fDz;
    fOuterOffset    = fRmax2 - fOuterSlope * fDz;
    fInnerTolerance = kConeTolerance * fSecRMin;
    fOuterTolerance = kConeTolerance * fSecRMax;

    if (fRmin2 > fRmin1) {
      fInnerConeApex     = 2 * fDz * fRmin1 / (fRmin2 - fRmin1);
      fTanInnerApexAngle = fRmin2 / (2 * fDz + fInnerConeApex);
    } else {
      fInnerConeApex     = 2 * fDz * fRmin2 / NonZero(fRmin1 - fRmin2);
      fTanInnerApexAngle = fRmin1 / (2 * fDz + fInnerConeApex);
    }
    if (fRmin1 == 0. || fRmin2 == 0.) fInnerConeApex = 0.;
    if (fRmin1 == 0.) fTanInnerApexAngle = fRmin2 / (2 * fDz);
    if (fRmin2 == 0.) fTanInnerApexAngle = fRmin1 / (2 * fDz);

    if (fRmax2 > fRmax1) {
      fOuterConeApex     = 2 * fDz * fRmax1 / (fRmax2 - fRmax1);
      fTanOuterApexAngle = fRmax2 / (2 * fDz + fOuterConeApex);
    } else {
      fOuterConeApex     = 2 * fDz * fRmax2 / NonZero(fRmax1 - fRmax2);
      fTanOuterApexAngle = fRmax1 / (2 * fDz + fOuterConeApex);
    }
    if (fRmax1 == 0. || fRmax2 == 0.) fOuterConeApex = 0.;
    if (fRmax1 == 0.) fTanOuterApexAngle = fRmax2 / (2 * fDz);
    if (fRmax2 == 0.) fTanOuterApexAngle = fRmax1 / (2 * fDz);

    fZNormInner = fTanRMin / NonZero(fSecRMin);
    fZNormOuter = -fTanRMax / NonZero(fSecRMax);
    fTolIz      = fDz - kHalfTolerance;
    fTolOz      = fDz + kHalfTolerance;
  }

  VECCORE_ATT_HOST_DEVICE
  void Init(T rmin1, T rmax1, T rmin2, T rmax2, T dz)
  {
    fRmin1 = rmin1 < 0. ? 0. : rmin1;
    fRmax1 = rmax1;
    fRmin2 = rmin2 < 0. ? 0. : rmin2;
    fRmax2 = rmax2;
    fDz    = dz;
    CalculateCached();
  }

  /// @brief Accumulate the outward normals of all incident radial and z surfaces.
  /// @details The returned vector is an unnormalized sum so callers can append
  /// normals from surfaces represented outside this radial section.
  /// @param p Point in the section's local coordinates.
  /// @param[out] norm Raw incident-surface sum, or the nearest-surface direction
  /// when the point is not on a surface.
  /// @param[out] count Number of incident surfaces included in `norm`; zero when
  /// the nearest-surface fallback is returned.
  /// @return True for one or two incident surfaces; false for no surface or more
  /// than two incident surfaces.
  VECCORE_ATT_HOST_DEVICE
  bool NormalWithCount(Vector3D<Precision> const &p, Vector3D<Precision> &norm, int &count) const
  {
    const Precision rho = Sqrt(p.Perp2());
    Precision distZ     = Abs(Abs(p.z()) - fDz);
    Precision distRMin;
    Precision distRMax;

    const Precision pRMin   = rho - p.z() * fTanRMin;
    const Precision widRMin = fRmin2 - fDz * fTanRMin;
    if (Abs(fRmin1 - fRmin2) < fInnerTolerance)
      distRMin = Abs(rho - fRmin2);
    else
      distRMin = Abs((pRMin - widRMin) / fSecRMin);

    const Precision pRMax   = rho - p.z() * fTanRMax;
    const Precision widRMax = fRmax2 - fDz * fTanRMax;
    if (Abs(fRmax1 - fRmax2) < fOuterTolerance)
      distRMax = Abs(rho - fRmax2);
    else
      distRMax = Abs((pRMax - widRMax) / fSecRMax);

    Vector3D<Precision> radialOuter(0., 0., fZNormOuter);
    Vector3D<Precision> radialInner(0., 0., fZNormInner);
    if (rho > kHalfTolerance) {
      radialOuter.Set(p.x() / rho / fSecRMax, p.y() / rho / fSecRMax, -fTanRMax / fSecRMax);
      radialInner.Set(-p.x() / rho / fSecRMin, -p.y() / rho / fSecRMin, fTanRMin / fSecRMin);
    }

    const bool radiallyInsideOuter = rho <= GetRadiusAtZ<false>(p.z()) + fOuterTolerance;
    const bool radiallyInsideInner =
        (fRmin1 == 0. && fRmin2 == 0.) || rho >= GetRadiusAtZ<true>(p.z()) - fInnerTolerance;
    const bool insideRadial = radiallyInsideOuter && radiallyInsideInner;

    count = 0;
    Vector3D<Precision> sum(0., 0., 0.);
    if (insideRadial && distZ <= kHalfTolerance) {
      const bool upper      = p.z() >= 0.;
      const bool degenerate = upper ? Abs(fRmax2 - fRmin2) < kTolerance : Abs(fRmax1 - fRmin1) < kTolerance;
      if (!degenerate) {
        sum += Vector3D<Precision>(0., 0., upper ? 1. : -1.);
        ++count;
      }
    }
    if (Abs(p.z()) <= fDz + kTolerance && distRMax <= fOuterTolerance) {
      sum += radialOuter;
      ++count;
    }
    if ((fRmin1 || fRmin2) && Abs(p.z()) <= fDz + kTolerance && distRMin <= fInnerTolerance) {
      sum += radialInner;
      ++count;
    }

    if (count) {
      norm = sum;
      return count <= 2;
    }
    norm = distZ < distRMax ? Vector3D<Precision>(0., 0., p.z() < 0. ? -1. : 1.) : radialOuter;
    if ((fRmin1 || fRmin2) && distRMin < Min(distZ, distRMax)) norm = radialInner;
    return false;
  }

  VECCORE_ATT_HOST_DEVICE
  bool Normal(Vector3D<Precision> const &p, Vector3D<Precision> &norm) const
  {
    int count        = 0;
    const bool valid = NormalWithCount(p, norm, count);
    if (count > 1) norm = norm.Unit();
    return valid;
  }

  template <bool Inner>
  VECCORE_ATT_HOST_DEVICE T GetRadiusAtZ(T z) const
  {
    return Inner ? T(fInnerSlope * z + fInnerOffset) : T(fOuterSlope * z + fOuterOffset);
  }

  VECCORE_ATT_HOST_DEVICE
  void Print() const
  {
    printf("RadialConeStruct : {rmin1 %.2f, rmax1 %.2f, rmin2 %.2f, rmax2 %.2f, dz %.2f}", fRmin1, fRmax1, fRmin2,
           fRmax2, fDz);
  }
};

} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif
