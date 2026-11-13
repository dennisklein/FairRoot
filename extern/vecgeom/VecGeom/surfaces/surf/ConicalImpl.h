#ifndef VECGEOM_CONICAL_IMPL_H
#define VECGEOM_CONICAL_IMPL_H

#include "VecGeom/surfaces/surf/SurfaceHelper.h"
#include "VecGeom/surfaces/base/Equations.h"

namespace vgbrep {

template <typename Real_t>
struct SurfaceHelper<SurfaceType::kConical, Real_t> {

  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  SurfaceHelper() {}

  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  /// @brief Inside half-space function
  /// @param point Point in local surface coordinates
  /// @param tol tolerance for determining if a point is inside or not
  /// @return True if the point is behind the normal within kTolerance (surface is included)
  bool Inside(Vector3D<Real_t> const &point, bool flip, Real_t radius, Real_t slope,
              Real_t tol = vecgeom::kToleranceStrict<Real_t>)
  {
    int flipsign      = (radius < 0) ? -1 : 1;
    int bool_flipsign = !flip ? 1 : -1;
    Real_t coneRadius = vecgeom::Abs(radius) + point.z() * slope;
    Real_t rho        = point.Perp();
    return flipsign * (rho - coneRadius) < bool_flipsign * tol;
  }

  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  /// @brief Find signed distance to next intersection from local point.
  /// @param point Point in local surface coordinates
  /// @param dir Direction in the local surface coordinates
  /// @param left_side Flag specifying if the surface is intersected from the left-side that defines the normal
  /// @param distance Computed distance to surface
  /// @return Validity of the intersection
  bool Intersect(Vector3D<Real_t> const &point, Vector3D<Real_t> const &dir, bool left_side, Real_t &distance,
                 bool &two_solutions, Real_t &safety, Real_t radius, Real_t slope)
  {
    QuadraticCoef<Real_t> coef;
    Real_t roots[2];
    int numroots      = 0;
    bool flip_exiting = left_side ^ (radius < 0);
    ConeEq<Real_t>(point, dir, vecgeom::Abs(radius), slope, coef);
    QuadraticSolver(coef, roots, numroots);
    two_solutions = (numroots == 2 && roots[0] > -vecgeom::kToleranceStrict<Real_t> &&
                     roots[1] > -vecgeom::kToleranceStrict<Real_t>);
    for (auto i = 0; i < numroots; ++i) {
      distance                = roots[i];
      Vector3D<Real_t> onsurf = point + distance * dir;
      // Exclude solutions beyond the tip of the cone.
      if (vecgeom::Abs(radius) + onsurf[2] * slope < 0) continue;
      Real_t rhosurf = vecgeom::Sqrt(onsurf[0] * onsurf[0] + onsurf[1] * onsurf[1]);
      bool entering;
      if (rhosurf > vecgeom::kToleranceStrict<Real_t>) {
        Vector3D<Real_t> normal(onsurf[0], onsurf[1], -rhosurf * slope);
        entering = dir.Dot(normal) < 0;
      } else {
        // The intersection is the tip itself, where the normal vanishes. The ray
        // enters iff it stays within the cone's half-angle on the side the cone
        // opens to, which the normal test cannot express.
        entering = slope * dir[2] > vecgeom::Sqrt(dir[0] * dir[0] + dir[1] * dir[1]);
      }
      bool hit = flip_exiting ^ entering;
      // First solution giving a valid hit wins
      if (hit) {
        if (distance < -vecgeom::kToleranceStrict<Real_t>) {
          Real_t rho     = point.Perp();
          auto distanceR = vecgeom::Abs(radius) + point[2] * slope - rho;
          Real_t calf    = static_cast<Real_t>(1) / vecgeom::Sqrt(static_cast<Real_t>(1) + slope * slope);
          safety         = distanceR * calf;
        }
        return true;
      }
    }
    // A ray through the tip meets the cone in a double root, which cancellation
    // in the discriminant can lose entirely, so look for that crossing directly.
    return IntersectTip(point, dir, flip_exiting, distance, radius, slope);
  }

  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  /// @brief Crossing of the cone tip, where the surface normal is undefined.
  /// @return True if the ray passes through the tip and continues on the solid side
  bool IntersectTip(Vector3D<Real_t> const &point, Vector3D<Real_t> const &dir, bool flip_exiting, Real_t &distance,
                    Real_t radius, Real_t slope)
  {
    if (slope == Real_t(0)) return false;
    // tip of the cone on the axis, where its radius vanishes
    Vector3D<Real_t> tip(0, 0, -vecgeom::Abs(radius) / slope);
    Vector3D<Real_t> totip = tip - point;
    Real_t dtip            = totip.Dot(dir);
    if ((totip - dtip * dir).Mag() > vecgeom::kToleranceStrict<Real_t>) return false;
    bool entering = slope * dir[2] > vecgeom::Sqrt(dir[0] * dir[0] + dir[1] * dir[1]);
    if (!(flip_exiting ^ entering)) return false;
    distance = dtip;
    return true;
  }

  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  /// @brief Computes the isotropic safe distance to unplaced surfaces
  /// @param point Point in local surface coordinates
  /// @param left_side Flag specifying if the surface is intersected from the left-side that defines the normal
  /// @param distance Computed isotropic safety
  /// @param onsurf Projection of the point on surface
  /// @return Validity of the calculation
  bool Safety(Vector3D<Real_t> const &point, bool left_side, Real_t &distance, Vector3D<Real_t> &onsurf, Real_t radius,
              Real_t slope) const
  {
    Real_t coneRadius = vecgeom::Abs(radius) + point[2] * slope;
    Real_t rho        = point.Perp();
    bool flip_exiting = left_side ^ (radius < 0);
    auto distanceR    = flip_exiting ? coneRadius - rho : rho - coneRadius;
    Real_t calf       = static_cast<Real_t>(1) / vecgeom::Sqrt(static_cast<Real_t>(1) + slope * slope);
    distance          = distanceR * calf;
    // We only use for the ZPhi frame safety the z of the point propagated on the cone surface
    onsurf = point;
    onsurf[2] += distance * slope * calf;
    return distance > -vecgeom::kToleranceDist<Real_t>;
  }
};

} // namespace vgbrep

#endif
