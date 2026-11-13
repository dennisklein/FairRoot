/*
 * ConeStruct.h
 *
 *  Created on: May 09, 2017
 *      Author: Raman Sehgal
 */
#ifndef VECGEOM_CONESTRUCT_H_
#define VECGEOM_CONESTRUCT_H_

#include "VecGeom/base/Global.h"
#include "VecGeom/volumes/Wedge_Evolution.h"
#include "VecGeom/management/Logger.h"

namespace vecgeom {

VECGEOM_DEVICE_DECLARE_CONV_TEMPLATE(struct, ConeStruct, typename);

inline namespace VECGEOM_IMPL_NAMESPACE {

// a plain and lightweight struct to encapsulate data members of a Cone
template <typename T = double>
struct ConeStruct {
  // Cone defining parameters
  T fRmin1{0.};
  T fRmax1{0.};
  T fRmin2{0.};
  T fRmax2{0.};
  T fDz{0.};
  T fSPhi{0.};
  T fDPhi{0.};

  evolution::Wedge fPhiWedge;

  // vectors characterizing the normals of phi planes
  // makes task to detect phi sektors very efficient
  Vector3D<Precision> fNormalPhi1;
  Vector3D<Precision> fNormalPhi2;
  Precision fAlongPhi1x{0.};
  Precision fAlongPhi1y{0.};
  Precision fAlongPhi2x{0.};
  Precision fAlongPhi2y{0.};

  // Some Cached value, try to reduce them
  // Some precomputed values to avoid divisions etc
  Precision fInnerSlope{0.}; // "gradient" of inner surface in z direction
  Precision fOuterSlope{0.}; // "gradient" of outer surface in z direction
  Precision fInnerOffset{0.};
  Precision fOuterOffset{0.};
  Precision fInnerTolerance{0.}; // tolerance on radial direction for inner surface
  Precision fOuterTolerance{0.}; // tolerance on radial direction for outer surface
  // Values to be cached
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
  Precision Capacity() const
  {
    return (fDz * fDPhi / 3.) *
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

    // check this very carefully
    fInnerSlope     = -(fRmin1 - fRmin2) / (2. * fDz);
    fOuterSlope     = -(fRmax1 - fRmax2) / (2. * fDz);
    fInnerOffset    = fRmin2 - fInnerSlope * fDz;
    fOuterOffset    = fRmax2 - fOuterSlope * fDz;
    fInnerTolerance = kConeTolerance * fSecRMin;
    fOuterTolerance = kConeTolerance * fSecRMax;

    if (fRmin2 > fRmin1) {
      fInnerConeApex     = 2 * fDz * fRmin1 / (fRmin2 - fRmin1);
      fTanInnerApexAngle = fRmin2 / (2 * fDz + fInnerConeApex);
    } else { // Should we add a check if(fRmin1 > fRmin2)
      fInnerConeApex     = 2 * fDz * fRmin2 / NonZero(fRmin1 - fRmin2);
      fTanInnerApexAngle = fRmin1 / (2 * fDz + fInnerConeApex);
    }

    if (fRmin1 == 0. || fRmin2 == 0.) fInnerConeApex = 0.;

    if (fRmin1 == 0.) fTanInnerApexAngle = fRmin2 / (2 * fDz);
    if (fRmin2 == 0.) fTanInnerApexAngle = fRmin1 / (2 * fDz);

    if (fRmax2 > fRmax1) {
      fOuterConeApex     = 2 * fDz * fRmax1 / (fRmax2 - fRmax1);
      fTanOuterApexAngle = fRmax2 / (2 * fDz + fOuterConeApex);
    } else { // Should we add a check if(fRmax1 > fRmax2)
      fOuterConeApex     = 2 * fDz * fRmax2 / NonZero(fRmax1 - fRmax2);
      fTanOuterApexAngle = fRmax1 / (2 * fDz + fOuterConeApex);
    }

    if (fRmax1 == 0. || fRmax2 == 0.) fOuterConeApex = 0.;

    if (fRmax1 == 0.) fTanOuterApexAngle = fRmax2 / (2 * fDz);
    if (fRmax2 == 0.) fTanOuterApexAngle = fRmax1 / (2 * fDz);

    fZNormInner = fTanRMin / NonZero(fSecRMin);
    fZNormOuter = -fTanRMax / NonZero(fSecRMax);

    fTolIz = fDz - kHalfTolerance;
    fTolOz = fDz + kHalfTolerance;

    // DetectConvexity();
  }

  VECCORE_ATT_HOST_DEVICE
  void Print() const
  {
    printf("ConeStruct :  {rmin1 %.2f, rmax1 %.2f, rmin2 %.2f, "
           "rmax2 %.2f, dz %.2f, phistart %.2f, deltaphi %.2f}",
           fRmin1, fRmax1, fRmin2, fRmax2, fDz, fSPhi, fDPhi);
  }

  void Print(std::ostream &os) const { os << "UnplacedCone; please implement Print to outstream\n"; }

  VECCORE_ATT_HOST_DEVICE
  bool IsFullPhi() const { return fDPhi == kTwoPi; }

  VECCORE_ATT_HOST_DEVICE
  bool Normal(Vector3D<Precision> const &p, Vector3D<Precision> &norm) const
  {
    int noSurfaces = 0;
    Precision rho, pPhi;
    Precision distZ, distRMin, distRMax;
    Precision distSPhi = kInfLength, distEPhi = kInfLength;
    Precision pRMin, widRMin;
    Precision pRMax, widRMax;

    // const double kHalfTolerance = 0.5 * kTolerance;

    Vector3D<Precision> sumnorm(0., 0., 0.), nZ = Vector3D<Precision>(0., 0., 1.);
    Vector3D<Precision> nR, nr(0., 0., 0.), nPs, nPe;
    norm = sumnorm;

    // do not use an extra fabs here -- negative/positive distZ tells us when point is outside or inside
    distZ = Abs(p.z()) - fDz;
    rho   = Sqrt(p.x() * p.x() + p.y() * p.y());

    pRMin   = rho - p.z() * fTanRMin;
    widRMin = fRmin2 - fDz * fTanRMin;
    if (Abs(fRmin1 - fRmin2) < fInnerTolerance)
      distRMin = (rho - fRmin2);
    else
      distRMin = (pRMin - widRMin) / fSecRMin;

    pRMax   = rho - p.z() * fTanRMax;
    widRMax = fRmax2 - fDz * fTanRMax;
    if (Abs(fRmax1 - fRmax2) < fOuterTolerance)
      distRMax = (rho - fRmax2);
    else
      distRMax = (pRMax - widRMax) / fSecRMax;

    bool inside = false;
    if (distZ < kTolerance) inside = distRMax < fOuterTolerance;
    if (inside && (fRmin1 || fRmin2)) inside = distRMin > -fInnerTolerance;

    distZ    = Abs(distZ);
    distRMax = Abs(distRMax);
    distRMin = Abs(distRMin);

    // keep track of nearest normal, needed in case point is not on a surface
    Precision distNearest           = distZ;
    Vector3D<Precision> normNearest = nZ;
    if (p.z() < 0.) normNearest.Set(0, 0, -1.);

    if (!IsFullPhi()) {
      if (rho) { // Protected against (0,0,z)
        pPhi = ATan2(p.y(), p.x());

        if (pPhi < fSPhi - kHalfTolerance)
          pPhi += 2 * kPi;
        else if (pPhi > fSPhi + fDPhi + kHalfTolerance)
          pPhi -= 2 * kPi;

        distSPhi = rho * (pPhi - fSPhi);
        distEPhi = rho * (pPhi - fSPhi - fDPhi);
        inside   = inside && (distSPhi > -kTolerance) && (distEPhi < kTolerance);
        distSPhi = Abs(distSPhi);
        distEPhi = Abs(distEPhi);
      }

      else if (!(fRmin1) || !(fRmin2)) {
        distSPhi = 0.;
        distEPhi = 0.;
      }
      nPs = Vector3D<Precision>(fAlongPhi1y, -fAlongPhi1x, 0);
      nPe = Vector3D<Precision>(-fAlongPhi2y, fAlongPhi2x, 0);
    }

    if (rho > kHalfTolerance) {
      // The cone slopes are fixed; only the point's radial direction needs normalization.
      const Precision invRho     = 1. / rho;
      const Precision outerScale = invRho * fInvSecRMax;
      nR                         = Vector3D<Precision>(p.x() * outerScale, p.y() * outerScale, fZNormOuter);
      if (fRmin1 || fRmin2) {
        const Precision innerScale = invRho * fInvSecRMin;
        nr                         = Vector3D<Precision>(-p.x() * innerScale, -p.y() * innerScale, fZNormInner);
      }
    }

    if (inside && distZ <= kHalfTolerance) {
      const bool onUpperEnd = p.z() >= 0.;
      const bool degenerateEndRing =
          onUpperEnd ? (Abs(fRmax2 - fRmin2) < kTolerance) : (Abs(fRmax1 - fRmin1) < kTolerance);
      if (!degenerateEndRing) {
        noSurfaces++;
        if (onUpperEnd)
          sumnorm += nZ;
        else
          sumnorm.Set(0, 0, -1.);
      }
    }

    if (inside && distRMax <= fOuterTolerance) {
      noSurfaces++;
      sumnorm += nR;
    } else if (noSurfaces == 0 && distRMax < distNearest) {
      distNearest = distRMax;
      normNearest = nR;
    }

    if (fRmin1 || fRmin2) {
      if (inside && distRMin <= fInnerTolerance) {
        noSurfaces++;
        sumnorm += nr;
      } else if (noSurfaces == 0 && distRMin < distNearest) {
        distNearest = distRMin;
        normNearest = nr;
      }
    }

    if (!IsFullPhi()) {
      if (inside && distSPhi <= kHalfTolerance) {
        noSurfaces++;
        sumnorm += nPs;
      } else if (noSurfaces == 0 && distSPhi < distNearest) {
        distNearest = distSPhi;
        normNearest = nPs;
      }
      if (inside && distEPhi <= kHalfTolerance) {
        noSurfaces++;
        sumnorm += nPe;
      } else if (noSurfaces == 0 && distEPhi < distNearest) {
        // No more check on distNearest, no need to assign to it.
        // distNearest = distEPhi;
        normNearest = nPe;
      }
    }
    // Final checks
    if (noSurfaces == 0)
      norm = normNearest;
    else if (noSurfaces == 1)
      norm = sumnorm;
    else
      norm = sumnorm.Unit();

    bool valid = noSurfaces != 0;
    if (noSurfaces > 2) {
      // return valid=false for noSurfaces > 2
      valid = false;
    }

    return valid;
  }

  VECCORE_ATT_HOST_DEVICE
  void SetAndCheckSPhiAngle(Precision sPhi)
  {
    // Full phi has no meaningful start angle.
    if (fDPhi == kTwoPi) {
      fSPhi = 0;
    } else {
      // Ensure fSphi in 0-2PI or -2PI-0 range if shape crosses 0
      if (sPhi < 0) {
        fSPhi = kTwoPi - Fmod(Abs(sPhi), kTwoPi);
      } else {
        fSPhi = Fmod(sPhi, kTwoPi);
      }
      if (fSPhi + fDPhi > kTwoPi) {
        fSPhi -= kTwoPi;
      }
    }

    // Update Wedge
    fPhiWedge.Set(fDPhi, fSPhi);
    fPhiWedge.UpdateNormals();
    // Update cached values.
    GetAlongVectorToPhiSector(fSPhi, fAlongPhi1x, fAlongPhi1y);
    GetAlongVectorToPhiSector(fSPhi + fDPhi, fAlongPhi2x, fAlongPhi2y);
  }

  VECCORE_ATT_HOST_DEVICE
  void SetAndCheckDPhiAngle(Precision dPhi)
  {
    if (dPhi <= 0 || dPhi >= kTwoPi - 0.5 * kAngTolerance) {
      fDPhi = kTwoPi;
      fSPhi = 0;
    } else {
      fDPhi = dPhi;
    }
    // Update Wedge
    fPhiWedge.Set(fDPhi, fSPhi);
    fPhiWedge.UpdateNormals();
    // Update cached values.
    GetAlongVectorToPhiSector(fSPhi, fAlongPhi1x, fAlongPhi1y);
    GetAlongVectorToPhiSector(fSPhi + fDPhi, fAlongPhi2x, fAlongPhi2y);
  }

  VECCORE_ATT_HOST_DEVICE
  static void GetAlongVectorToPhiSector(Precision phi, Precision &x, Precision &y)
  {
    x = Cos(phi);
    y = Sin(phi);
  }

  void SetRmin1(Precision const &arg)
  {
    fRmin1 = arg;
    CalculateCached();
  }
  void SetRmax1(Precision const &arg)
  {
    fRmax1 = arg;
    CalculateCached();
  }
  void SetRmin2(Precision const &arg)
  {
    fRmin2 = arg;
    CalculateCached();
  }
  void SetRmax2(Precision const &arg)
  {
    fRmax2 = arg;
    CalculateCached();
  }
  void SetDz(Precision const &arg)
  {
    fDz = arg;
    CalculateCached();
  }
  void SetSPhi(Precision const &arg)
  {
    fSPhi = arg;
    SetAndCheckSPhiAngle(fSPhi);
    // DetectConvexity();
  }
  void SetDPhi(Precision const &arg)
  {
    fDPhi = arg;
    SetAndCheckDPhiAngle(fDPhi);
    // DetectConvexity();
  }

  VECCORE_ATT_HOST_DEVICE
  Precision GetTolIz() const { return fTolIz; }
  VECCORE_ATT_HOST_DEVICE
  Precision GetTolOz() const { return fTolOz; }

  VECCORE_ATT_HOST_DEVICE
  evolution::Wedge const &GetWedge() const { return fPhiWedge; }

  // constructors
  ConeStruct() = default;

  VECCORE_ATT_HOST_DEVICE
  ConeStruct(T const &rmin1, T const &rmax1, T const &rmin2, T const &rmax2, T const &z, T const &sphi, T const &dphi)
  {
    Init(rmin1, rmax1, rmin2, rmax2, z, sphi, dphi);
  }

  VECCORE_ATT_HOST_DEVICE
  void Init(T const &rmin1, T const &rmax1, T const &rmin2, T const &rmax2, T const &z, T const &sphi, T const &dphi)
  {
    fRmin1 = rmin1 < 0.0 ? 0.0 : rmin1;
    fRmax1 = rmax1;
    fRmin2 = rmin2 < 0.0 ? 0.0 : rmin2;
    fRmax2 = rmax2;
    fDz    = z;
    fSPhi  = sphi;
    fDPhi  = dphi;
    fPhiWedge.Init(dphi, sphi);
    SetAndCheckDPhiAngle(dphi);
    SetAndCheckSPhiAngle(sphi);
    CalculateCached();
    // DetectConvexity();
  }
};
} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif
