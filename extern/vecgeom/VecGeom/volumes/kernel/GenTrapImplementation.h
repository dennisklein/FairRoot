/*
 * \file GenTrapImplementation.h
 * \brief Navigation kernels for the Generic Trapezoid (Arb8).
 *
 * This header implements the runtime kernels (Contains/Inside, Distance, Safety, Normal)
 * that operate on the immutable data stored in GenTrapStruct (see GenTrapStruct.h).
 *
 * Performance notes:
 *  - Planar cases can be short-circuited via TessellatedSection.
 *  - Lateral twisted faces (bilinear patches) use a cheap v(z)-window
 *    to cull impossible distance ranges before doing any quadratic/UV work.
 *  - For ray/face intersection we adopt a half-open step window [lmin, lmax),
 *    so hits at exactly lmax are considered "not in range" and lmin is inclusive
 *    up to a small negative tolerance to absorb FP noise.
 *
 *  Created on: Aug 2, 2014
 *      Author: swenzel
 *   Review/completion: Nov 4, 2015
 *      Author: mgheata
 *  Full scalar refactoring: andrei.gheata@cern.ch
 */

#ifndef VECGEOM_VOLUMES_KERNEL_GENTRAPIMPLEMENTATION_H_
#define VECGEOM_VOLUMES_KERNEL_GENTRAPIMPLEMENTATION_H_

#include "VecGeom/base/Global.h"
#include "VecGeom/base/Vector2D.h"

#include "VecGeom/volumes/kernel/GenericKernels.h"
#include "VecGeom/volumes/kernel/BoxImplementation.h"
#include "VecGeom/volumes/SurfaceHitView.h"
#include "VecGeom/volumes/UnplacedGenTrap.h"

#include <iostream>

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE(struct GenTrapImplementation;);
VECGEOM_DEVICE_DECLARE_CONV(struct, GenTrapImplementation);

inline namespace VECGEOM_IMPL_NAMESPACE {

class PlacedGenTrap;
class UnplacedGenTrap;

template <typename T>
struct GenTrapStruct;

struct GenTrapImplementation {

  using Vertex_t         = Vector3D<Precision>;
  using PlacedShape_t    = PlacedGenTrap;
  using UnplacedStruct_t = GenTrapStruct<Precision>;
  using UnplacedVolume_t = UnplacedGenTrap;

  /// Boundary selected by the existing DistanceToOut minimization.
  enum class ExitSurface : signed char { kNone = -1, kZMin, kZMax, kLateral0, kLateral1, kLateral2, kLateral3 };

  /**
   * \brief Point containment test (fast path, no surface classification).
   *
   * Checks whether a point lies inside the generalized trapezoid volume.
   * This is typically cheaper than Inside() as it does not distinguish
   * between kInside and kSurface -- it only returns a boolean.
   *
   * \tparam Real_v   Scalar floating type
   * \param unplaced  Geometry data structure of the GenTrap
   * \param point     Query point in local coordinates
   * \param[out] inside  True if the point is not rejected by any face
   */
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void Contains(UnplacedStruct_t const &unplaced,
                                                                    Vector3D<Real_v> const &point, bool &inside);

  /**
   * \brief Classify a point with respect to the solid.
   *
   * Returns one of EInside::{kInside, kSurface, kOutside}. This may perform
   * additional tolerance checks vs. Contains().
   *
   * \tparam Real_v     Scalar or SIMD floating type
   * \param unplaced    Geometry data structure of the GenTrap
   * \param point       Query point in local coordinates
   * \param[out] inside Classification result
   */
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void Inside(UnplacedStruct_t const &unplaced,
                                                                  Vector3D<Real_v> const &point, Inside_t &inside);

  /**
   * \brief Shared kernel used by Contains() and Inside().
   *
   * Evaluates the bilinear/planar surface equations on the four lateral faces
   * (and the implicit Z slabs) to determine quick-reject/accept masks.
   * When ForInside=true, it additionally computes a conservative
   * \c completelyinside mask.
   *
   * \tparam Real_v   Scalar type (assumed scalar in this implementation)
   * \tparam ForInside If true, also compute \c completelyinside; otherwise only outside mask
   * \param unplaced   Geometry data structure of the GenTrap
   * \param point      Query point
   * \param[out] completelyinside True if strictly inside all faces by tolerance
   * \param[out] completelyoutside True if outside any face by tolerance
   */
  template <typename Real_v, bool ForInside>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void GenericKernelForContainsAndInside(
      UnplacedStruct_t const &unplaced, Vector3D<Real_v> const &point, bool &completelyinside, bool &completelyoutside);

  /**
   * \brief Ray entry distance from an external point.
   *
   * Computes the distance along \p direction from \p point to enter the solid,
   * honoring \p stepMax. Returns \c -1 for points that are already inside.
   *
   * \tparam Real_v Scalar or SIMD floating type
   * \param unplaced Geometry data structure of the GenTrap
   * \param point    Start point (expected outside)
   * \param direction Normalized or non-normalized direction; sign matters
   * \param stepMax  Upper bound for the step; intersections beyond are ignored
   * \param[out] distance The entry distance, or +inf if no hit, or -1 if already inside
   */
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void DistanceToIn(UnplacedStruct_t const &unplaced,
                                                                        Vector3D<Real_v> const &point,
                                                                        Vector3D<Real_v> const &direction,
                                                                        Real_v const &stepMax, Real_v &distance);

  /**
   * \brief Ray exit distance with optional supporting-normal output.
   *
   * Computes the distance from \p point along \p direction to leave the solid.
   *
   * Caps and planar lateral faces support the whole solid. Twisted lateral
   * exits conservatively leave hitInfo->fLeftBehind false and normal storage
   * untouched.
   * A null hitInfo retains the ordinary distance-only path.
   *
   * \tparam Real_v Scalar floating-point type
   * \param unplaced Geometry data structure of the GenTrap
   * \param point Start point inside or on the surface
   * \param direction Unit ray direction
   * \param stepMax  Currently unused upper bound (kept for API symmetry)
   * \param[out] distance Exit distance; tolerated surface starts may yield a
   * small negative value. Outside starts return -1, or -kTolerance in the planar helper.
   * \param[out] hitInfo Optional output; cleared first, with a unit normal
   * written only when requested and the selected boundary supports the solid
   */
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void DistanceToOut(UnplacedStruct_t const &unplaced,
                                                                         Vector3D<Real_v> const &point,
                                                                         Vector3D<Real_v> const &direction,
                                                                         Real_v const &stepMax, Real_v &distance,
                                                                         SurfaceHitView<Real_v> *hitInfo = nullptr);

  /** @brief Compute the exit distance with compile-time optional winner recording. */
  template <typename Real_v, bool RecordHitInfo>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void DistanceToOutImpl(UnplacedStruct_t const &unplaced,
                                                                             Vector3D<Real_v> const &point,
                                                                             Vector3D<Real_v> const &direction,
                                                                             Real_v const &stepMax, Real_v &distance,
                                                                             SurfaceHitView<Real_v> *hitInfo);

  /** @brief Convert a selected cap or planar lateral face into exit-normal output. */
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SetExitHit(UnplacedStruct_t const &unplaced,
                                                                      ExitSurface surface,
                                                                      SurfaceHitView<Real_v> *hitInfo);

  /**
   * \brief Conservative safety from outside to the solid.
   *
   * Returns an underestimate of the distance one must travel to reach the surface
   * from \p point. If \p point is outside the bounding box, the box safety is used.
   *
   * \tparam Real_v Scalar type (assumed scalar here)
   * \param unplaced Geometry data structure of the GenTrap
   * \param point    Query point
   * \param[out] safety Non-negative lower bound to the closest surface
   */
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SafetyToIn(UnplacedStruct_t const &unplaced,
                                                                      Vector3D<Real_v> const &point, Real_v &safety);

  /**
   * \brief Conservative safety from inside to the exterior.
   *
   * Returns an underestimate to the nearest exit surface.
   * If the point is found outside along Z, returns -1.
   *
   * \tparam Real_v Scalar type (assumed scalar here)
   * \param unplaced Geometry data structure of the GenTrap
   * \param point    Query point
   * \param[out] safety Lower bound to the nearest exit; -1 if already outside on Z
   */
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void SafetyToOut(UnplacedStruct_t const &unplaced,
                                                                       Vector3D<Real_v> const &point, Real_v &safety);

  /**
   * \brief Compute a unit-length outward normal at a surface point.
   *
   * Caps take precedence within tolerance; otherwise incident lateral normals
   * are averaged. Unmatched points use the nearest-edge approximation, with
   * twisted-face validity determined by tolerance-extended patch coordinates.
   *
   * \tparam Real_v Scalar type (assumed scalar here)
   * \param unplaced Geometry data structure of the GenTrap
   * \param point    Point on (or near) a surface
   * \param[out] normal Unit outward normal (undefined if !valid)
   * \param[out] valid  True if a consistent normal could be determined
   */
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void NormalKernel(UnplacedStruct_t const &unplaced,
                                                                        Vector3D<Real_v> const &point,
                                                                        Vector3D<Real_v> &normal, bool &valid);

  /// Utility functions

  /**
   * \brief Distance computation to a lateral surface.
   *
   * Handles both planar and twisted bilinear faces. For planar faces it reduces
   * to a plane intersection followed by in-range checks; for twisted faces it
   * solves a quadratic derived from the bilinear patch intersection.
   *
   * \tparam Real_v Precision type
   * \param unplaced Geometry data structure of the GenTrap
   * \param point Point
   * \param dir Direction (need not be unit length)
   * \param i Surface index [0,3]
   * \param in Inside state of the point (true if the ray starts inside)
   * \param lmin Minimum distance limit used to early-reject intersections
   * \param lmax Maximum distance limit used to early-reject intersections
   * \param assumeInsideStart True when the caller guarantees an inside/surface exit start
   * \return Distance to surface (\c +inf if no accepted hit lies in the window)
   */
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static Real_v DistanceToSurf(
      UnplacedStruct_t const &unplaced, Vector3D<Real_v> const &point, Vector3D<Real_v> const &dir, int i, bool in,
      Real_v lmin = Real_v(0), Real_v lmax = InfinityLength<Real_v>(), bool assumeInsideStart = false);

  /** \brief Fill vertices of the section at a given Z.
   *
   * \tparam Real_v Precision type
   * \param unplaced Geometry data structure of the GenTrap
   * \param z Z coordinate at which to evaluate the section
   * \param[out] vertices Out-parameter with the 4 XY vertices of the section
   */
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void FillVerticesAtZ(UnplacedStruct_t const &unplaced, Real_v z,
                                                                           Vector2D<Real_v> vertices[4]);

  /**
   * \brief Return the index of the edge closest to a given point (projected in XY).
   *
   * \tparam Real_v Scalar type
   * \param point   Query point (only x,y used)
   * \param vertxy  The four section vertices in XY
   * \param[out] iseg Index of the closest edge in [0,3]
   */
  template <class Real_v>
  VECCORE_ATT_HOST_DEVICE static void GetClosestEdge(Vector3D<Real_v> const &point, Vector2D<Real_v> const vertxy[4],
                                                     int &iseg);

  /** \brief Point-in-quad test on a Z section.
   *
   * Evaluates if the XY of \p point lies inside the quadrilateral defined by
   * \p v on the plane z=const of the corresponding section.
   *
   * \tparam Real_v Scalar type
   * \param point Query point (z is ignored except for consistency)
   * \param v     The four vertices of the section at the same z
   * \return True if inside (within tolerance); false if degenerate or outside
   */
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static bool InsideSection(Vector3D<Real_v> const &point,
                                                                         Vector3D<Real_v> const v[4]);

  /**
   * \brief Check whether a point on lateral surface \p isurf is within its parametric bounds.
   *
   * Uses the inverse mapping (XYZ->(u,v)) and tests 0<=u,v<=1 with tolerance.
   *
   * \tparam Real_v Scalar type
   * \param unplaced Geometry data
   * \param point    Point assumed on the surface
   * \param isurf    Surface index [0,3]
   * \return True if the (u,v) parameters lie within [0,1]^2 (with tolerance)
   */
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static bool InSurfLimits(UnplacedStruct_t const &unplaced,
                                                                        Vector3D<Real_v> const &point, int isurf);

  /**
   * \brief Safety helper across the four lateral faces.
   *
   * Computes a conservative signed safety by aggregating the maximum of
   * signed distances (planar faces) or Lipschitz-continuous bounds (twisted faces).
   *
   * \tparam Real_v Scalar type
   * \param unplaced Geometry data
   * \param point    Query point
   * \param safmax   Initial safety bound (typically Z safety)
   * \param inside   If true, returns negative safety; otherwise positive
   * \return Signed conservative safety
   */
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static Real_v SafetyArb4(UnplacedStruct_t const &unplaced,
                                                                        Vector3D<Real_v> const &point, Real_v safmax,
                                                                        bool inside = true);

  /**
   * \brief Normal vector to a bilinear surface (un-normalized).
   *
   * The surface is parameterized as:
   * \f[
   * P(u,v) = (1-u)(1-v)P_1 + u(1-v)P_2 + (1-u)vP_3 + uvP_4
   * \f]
   * whose (non-normalized) normal is \f\vec n = \partial_u P \times \partial_v P\f.
   * The user may validate containment in the surface patch if (u,v) parameters lie within [0,1]^2 (with tolerance)
   *
   * \tparam Real_v Scalar type
   * \param unplaced Geometry data
   * \param point Point assumed on surface i
   * \param i Index of the lateral surface [0,3]
   * \param[out] unorm Un-normalized outward normal
   * \param[out] u shear coordinate at point.z()
   * \param[out] v generator coordinate
   */
  template <typename Real_v>
  VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE static void UNormal(UnplacedStruct_t const &unplaced,
                                                                   Vector3D<Real_v> const &point, int i,
                                                                   Vector3D<Real_v> &unorm, Real_v &u, Real_v &v);

  /** \brief Inverse mapping from cartesian point on a lateral face to (u,v).
   *
   * Implements the inverse of the bilinear patch equation restricted to the
   * generator pair (A,C) and (B,D) bounding the face.
   *
   * \tparam Real_v Scalar type
   * \param unplaced Geometry data
   * \param point Point in cartesian space, assumed on the bilinear surface i
   * \param i Index of the surface the point is on in the range [0, 3]
   * \return (u,v) coordinates within [0,1]^2 if the point is on the face
   */
  template <typename Real_v>
  VECCORE_ATT_HOST_DEVICE static Vector2D<Real_v> XYZtoUV(UnplacedStruct_t const &unplaced,
                                                          Vector3D<Real_v> const &point, int i);
}; // End struct GenTrapImplementation

//********************************
//**** implementations start here
//********************************/

//______________________________________________________________________________
// --- Cheap, safe interval culls for twisted faces ---
/**
 * \brief Compute the admissible t-interval such that the parametric v(t) stays in [0,1].
 *
 * The bilinear patch is parameterized with v determined by z:
 * \f[
 *   v(t) = \frac{z(t) + D_z}{2 D_z} = \frac{p_z + t\,d_z + D_z}{2 D_z}.
 * \f]
 * We pass \p Dz2 = 1/(2*D_z) to avoid divides. The result interval is intersected
 * with the user-provided forward limit \p limit and returned as [tmin,tmax].
 *
 * \param ptZ   start point z
 * \param dirZ  ray z-component
 * \param Dz2   precomputed 1/(2*Dz)
 * \param limit upper bound on t from the caller (e.g. stepMax)
 * \param[out] tmin lower bound of the feasible window
 * \param[out] tmax upper bound of the feasible window
 * \return false if the feasible window is empty
 */
template <typename Real_v>
VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE bool GenTrap_VInterval(Real_v ptZ, Real_v dirZ, Real_v Dz2, Real_v limit,
                                                                    Real_v &tmin, Real_v &tmax)
{
  // v(t) = (ptZ + Dz + dirZ*t) / (2*Dz) must stay in [0,1]
  const Real_v kv = dirZ * Dz2;              // slope of v(t)
  const Real_v v0 = ptZ * Dz2 + Real_v(0.5); // v at t=0
  Real_v t0 = Real_v(0), t1 = limit;         // user range
  if (Abs(kv) < kToleranceDist<Real_v>) {
    // v is (almost) constant along the ray
    if (v0 < Real_v(0) - kToleranceDist<Real_v> || v0 > Real_v(1) + kToleranceDist<Real_v>) return false;
    tmin = t0;
    tmax = t1;
    return true;
  }
  const Real_v t_at_0 = (Real_v(0) - v0) / kv;
  const Real_v t_at_1 = (Real_v(1) - v0) / kv;
  Real_v a = Min(t_at_0, t_at_1), b = Max(t_at_0, t_at_1);
  t0 = Max(t0, a);
  t1 = Min(t1, b);
  if (t1 < t0) return false;
  tmin = t0;
  tmax = t1;
  return true;
}

//______________________________________________________________________________
template <typename Real_v>
VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE Real_v GenTrapImplementation::DistanceToSurf(
    UnplacedStruct_t const &unplaced, Vector3D<Real_v> const &point, Vector3D<Real_v> const &dir, int i, bool in,
    Real_v lmin, Real_v lmax, bool assumeInsideStart)
{
  // Return a large sentinel for "no hit"
  Real_v big = InfinityLength<Real_v>();
  VECGEOM_ASSERT(!unplaced.IsDegenerated(i) && "DistanceToSurf called with degenerated face");

  // NOTE on window semantics:
  //   * We use a half-open window [lmin, lmax) to keep composition simple when
  //     iterating faces; a hit exactly at lmax is considered "out of range".
  //   * We allow a tiny negative slack at lmin to absorb FP error when starting
  //     very close to a surface: (t < lmin - tol) is rejected.
  //   * For inside->out queries, the first positive root is a valid exit
  //     as soon as it passes the window tests, hence early-return is safe.

  if (!unplaced.IsTwisted(i)) {
    // -------------------------
    // Planar face intersection
    // -------------------------

    // Signed plane safety: <0 means inside, >0 outside w.r.t. stored outward normal
    const Real_v safety = unplaced.fSurf[i].EvaluatePlanar(point); // n.p + g
    const Real_v dn     = unplaced.fSurf[i].DotPlaneNormal(dir);

    // Reject points on the wrong side or going in the wrong direction.
    //  - side_ok   : start is consistent with inside/outside state (tolerant)
    //  - facing_ok : direction crosses the plane in the correct sense
    // An exit must face outward: an inward near-parallel ray can have a large
    // backward root even from strictly inside, which must not be clamped to zero.
    const bool facing_ok = in ? dn > Real_v(0.) : dn < -kToleranceDist<Real_v>;
    const bool side_ok   = (in && safety < Real_v(kTolerance)) || (!in && safety > Real_v(-kTolerance));
    if (!facing_ok || !side_ok) return big;

    // Clamp tolerated backward crossings before the window check, so a surface
    // start is not discarded merely because its plane root lies behind it.
    Real_v snext = Max(-safety / NonZero(dn), Real_v(0.));
    // window culls before UV work
    if ((snext < lmin - kToleranceDist<Real_v>) || (snext >= lmax)) return big;

    // From strictly inside, the minimum over all faces and caps selects the
    // finite patch; an individual plane candidate need not recover its bounds.
    if (assumeInsideStart && safety < Real_v(-kTolerance)) return snext;

    // Validate the hit point lies within the finite bilinear bounds of the face
    return InSurfLimits<Real_v>(unplaced, point + snext * dir, i) ? snext : big;
  } else {
    // ---------------------------
    // Twisted bilinear intersection
    // ---------------------------
    // Intersection reduces to a quadratic in t; we test candidate roots against
    // the [lmin,lmax) window first, then validate against patch (u,v) bounds.
    // The code below uses directly the surface equation, but it is slower:
    // Real_v a, b, c;
    // unplaced.fSurf[i].RayQuadratic(point, dir, a, b, c);

#if GENTRAP_USE_HYBRID_COEFFS
    const auto pointXY = point.XY();
    const auto dirXY   = dir.XY();
    const auto z       = point[2];
    const auto dz      = dir[2];

    // Cached per-face differences
    const auto &dMid = unplaced.fDmid[i];
    const auto &dT   = unplaced.fDT[i];

    // Cached scalars (tiny memory footprint)
    const Real_v TT = unplaced.fTT[i];
    const Real_v MM = unplaced.fMM[i];
    const Real_v MT = unplaced.fMT[i];

    // Helper: z-component of 2D cross (kept inlined by the compiler)
    auto xz = [](auto const &a, auto const &b) -> Real_v { return a.x() * b.y() - a.y() * b.x(); };

    // Quadratic coefficients (hybrid precompute form)
    const Real_v dT_x_dir   = xz(dT, dirXY);
    const Real_v dT_x_point = xz(dT, pointXY);
    const Real_v dMid_x_dir = xz(dMid, dirXY);
    const Real_v dMid_x_pt  = xz(dMid, pointXY);

    // Share point-dependent terms in the quadratic's Horner form.
    const Real_v zTT        = z * TT;
    const Real_v pointTerms = dT_x_point + MT + zTT;
    Real_v a                = dz * (dT_x_dir + TT * dz);
    Real_v b                = dMid_x_dir + z * dT_x_dir + (pointTerms + zTT) * dz;
    Real_v c                = dMid_x_pt + MM + z * pointTerms;
#else
    int j = (i + 1) % 4;

    // Unpack a few frequently used quantities
    auto const pointXY = point.XY();  // project to XY
    auto const dirXY   = dir.XY();    // project direction to XY
    auto const t       = unplaced.fT; // shear vectors per vertex

    auto const s1 = unplaced.fMidXY[i] + t[i] * point[2]; // generator at z
    auto const s2 = unplaced.fMidXY[j] + t[j] * point[2];
    auto ds       = s2 - s1;

    // Quadratic coefficients (original on-the-fly derivation)
    Real_v a = ((t[j] - t[i]).CrossZ(dirXY) + t[i].CrossZ(t[j]) * dir[2]) * dir[2];
    Real_v b = ds.CrossZ(dirXY) + ((t[j] - t[i]).CrossZ(pointXY) + s1.CrossZ(t[j]) - s2.CrossZ(t[i])) * dir[2];
    Real_v c = ds.CrossZ(pointXY) + s1.CrossZ(s2);
#endif

    constexpr Real_v tolerance = kToleranceArb4<Real_v>;
    if (assumeInsideStart) {
      // For regular interior starts the first forward zero is this face's exit
      // candidate. Starts in the tolerated surface band additionally select an
      // outward crossing, so an inward ray skips its root at zero.
      const bool regularInterior = c < -tolerance;
      if (a == Real_v(0.)) {
        if (b == Real_v(0.)) return big;
        const Real_v root = -c / b;
        if (!regularInterior && b <= Real_v(0.)) return big;
        if (root < lmin - tolerance || root >= lmax) return big;
        const Real_v candidate = Max(root, Real_v(0.));
        return regularInterior || InSurfLimits(unplaced, point + candidate * dir, i) ? candidate : big;
      }

      const Real_v halfB = Real_v(0.5) * b;
      const Real_v delta = halfB * halfB - a * c;
      if (delta < Real_v(0.)) return big;
      const Real_v q = -halfB - CopySign(Sqrt(delta), halfB);
      if (regularInterior && delta > Real_v(0.) && lmin == Real_v(0.)) {
        // Negative q gives the first forward root c/q; otherwise q/a is
        // the only possible forward root. Keep repeated roots and shifted
        // windows on the general path below.
        const bool smallRoot = q < Real_v(0.);
        const Real_v root    = (smallRoot ? c : q) / (smallRoot ? NonZero(q) : a);
        return root >= Real_v(0.) && root < lmax ? root : big;
      }
      const Real_v r1       = q / a;
      const Real_v r2       = c / NonZero(q);
      Real_v root           = InfinityLength<Real_v>();
      const bool r1InWindow = r1 >= lmin || (r1 >= lmin - tolerance && b > Real_v(0.));
      const bool r2InWindow = r2 >= lmin || (r2 >= lmin - tolerance && b > Real_v(0.));
      if (r1InWindow && (regularInterior || Real_v(2.) * a * r1 + b > Real_v(0.))) root = r1;
      if (r2InWindow && (regularInterior || Real_v(2.) * a * r2 + b > Real_v(0.))) root = Min(root, r2);
      if (root >= lmax) return big;
      const Real_v candidate = Max(root, Real_v(0.));
      return regularInterior || InSurfLimits(unplaced, point + candidate * dir, i) ? candidate : big;
    }

    Real_v dist[2] = {-InfinityLength<Real_v>(), -InfinityLength<Real_v>()};
    auto nroots    = QuadraticSolver(a, b, c, dist);

    // The solver counts roots using the tighter generic tolerance but preserves
    // both. Recover a discarded root within the twisted-surface band only for
    // a non-grazing crossing of this finite face in the requested sense.
    const Real_v discardedRoot = nroots < 2 ? dist[1 - nroots] : -InfinityLength<Real_v>();
    if (discardedRoot > -tolerance && InSurfLimits<Real_v>(unplaced, point, i)) {
      const auto grad         = unplaced.fSurf[i].EvaluateGradient(point);
      const Real_v projection = dir.Dot(grad);
      const bool nonGrazing   = projection * projection > kToleranceSquared * grad.Mag2();
      if (nonGrazing && (in ^ (projection < Real_v(0.)))) return Real_v(0.);
    }

    // Prevent grazing rays when starting outside: two roots at ~0 within tolerance
    if (!in && nroots == 2 && Abs(dist[0]) < tolerance && Abs(dist[1]) < tolerance) return big;

    for (auto k = 2 - nroots; k < 2; ++k) {
      const Real_v troot = dist[k];

      // step window culls before UV/normal work
      if ((troot < lmin - kToleranceDist<Real_v>) || (troot >= lmax)) continue;

      // If the point is inside and not starting from a boundary, this is a minimization problem
      if (in && troot > tolerance) return troot; // inside: early minimum positive root is the exit

      // Allow roundoff at finite-patch edges, not a geometric tolerance in
      // dimensionless UV: the latter can admit remote hits on long faces.
      constexpr Real_v uvTolerance = Real_v(8.) * kEpsilonT<Real_v>;
      auto hit                     = point + troot * dir;
      auto uv                      = XYZtoUV(unplaced, hit, i);
      auto u                       = uv[0];
      auto v                       = uv[1];
      bool in_patch =
          (u > -uvTolerance) && (u < Real_v(1.) + uvTolerance) && (v > -uvTolerance) && (v < Real_v(1.) + uvTolerance);
      if (!in_patch) continue;
      // Crossing sense from implicit gradient (cheaper than building full geometric normal)
      auto grad       = unplaced.fSurf[i].EvaluateGradient(hit);
      bool crosses_ok = in ^ (dir.Dot(grad) < Real_v(0.));
      if (crosses_ok) return troot;
    }
  }
  return big;
}

//______________________________________________________________________________
template <typename Real_v, bool ForInside>
VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE void GenTrapImplementation::GenericKernelForContainsAndInside(
    UnplacedStruct_t const &unplaced, Vector3D<Real_v> const &point, bool &completelyinside, bool &completelyoutside)
{
  // First stage: AABB quick reject/accept; second stage: evaluate lateral faces.
  // For ForInside=true we also compute a strict "completelyinside" mask
  // by testing f_i(point) < -tolerance for all active faces. Points within
  // tolerance of the Z slabs are treated specially so we do not incorrectly
  // classify them as strictly inside (we then return kSurface in Inside()).

  // Quick reject/accept using the axis-aligned bounding box (AABB)
  Vector3D<Real_v> halfsize(unplaced.fBBdimensions[0], unplaced.fBBdimensions[1], unplaced.fBBdimensions[2]);
  BoxImplementation::GenericKernelForContainsAndInside<Real_v, true>(halfsize, point - unplaced.fBBorigin,
                                                                     completelyinside, completelyoutside);

  constexpr Real_v tolerance = kToleranceArb4<Real_v>;
  if (completelyoutside) return; // outside the AABB => outside the solid

  // If the point lies on the top/bottom faces within tolerance, skip tightening of the inside mask
  auto on_surfz = (unplaced.fDz - Abs(point.z())) < Real_v(kTolerance);

  // Evaluate implicit equations for lateral faces; any positive value means "outside that face"
  VECGEOM_PRAGMA_UNROLL(4)
  for (int i = 0; i < 4; i++) {
    if (unplaced.IsDegenerated(i)) continue;
    auto eqeval = unplaced.IsTwisted(i) ? unplaced.fSurf[i].Evaluate(point) : unplaced.fSurf[i].EvaluatePlanar(point);
    completelyoutside = eqeval > tolerance;
    if (completelyoutside) {
      completelyinside = false; // single face suffices to reject
      return;
    }
    if (ForInside && !on_surfz) completelyinside &= eqeval < -tolerance; // strictly inside all faces
  }
}

//______________________________________________________________________________
template <typename Real_v>
VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE void GenTrapImplementation::Contains(UnplacedStruct_t const &unplaced,
                                                                                  Vector3D<Real_v> const &point,
                                                                                  bool &inside)
{
  // Use tessellated helper when available (planar case optimization).
  // The helper is built on host and copied to device with scalar facet data.
  // This path handles all lateral faces as triangles/quads,
  // and is typically faster than evaluating bilinear patch equations
  // when the shape is purely planar. Z-slab rejection is kept explicit
  // to avoid tessellated work when clearly outside along Z.
  if (unplaced.fTslHelper) {
    // Check Z range
    inside = false;
    if (Abs(point.z()) > unplaced.fDz) return; // rejected by Z slab
    inside = unplaced.fTslHelper->Contains(point);
    return;
  }
  // Generic fallback
  bool unused  = false;
  bool outside = false;
  GenericKernelForContainsAndInside<Real_v, false>(unplaced, point, unused, outside);
  inside = !outside;
}

//______________________________________________________________________________
template <typename Real_v>
VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE void GenTrapImplementation::Inside(UnplacedStruct_t const &unplaced,
                                                                                Vector3D<Real_v> const &point,
                                                                                Inside_t &inside)
{
  if (unplaced.fTslHelper) {
    // Planar-case dispatcher with tessellated helper
    inside         = EInside::kOutside;
    Precision safZ = Abs(point.z()) - unplaced.fDz;
    if (safZ > kTolerance) return;     // definitely outside along Z
    bool insideZ = safZ < -kTolerance; // strictly between the Z planes
    inside       = unplaced.fTslHelper->Inside(point);
    if (insideZ || inside == EInside::kOutside) return; // kInside is already strict or outside
    inside = EInside::kSurface;                         // on the faces within tolerance
    return;
  }
  bool completelyinside;
  bool completelyoutside;
  GenericKernelForContainsAndInside<Real_v, true>(unplaced, point, completelyinside, completelyoutside);

  // Map masks to EInside classification
  inside = EInside::kSurface;
  if (completelyoutside) inside = EInside::kOutside;
  if (completelyinside) inside = EInside::kInside;
}

//______________________________________________________________________________
template <typename Real_v>
VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE void GenTrapImplementation::DistanceToIn(UnplacedStruct_t const &unplaced,
                                                                                      Vector3D<Real_v> const &point,
                                                                                      Vector3D<Real_v> const &direction,
                                                                                      Real_v const &stepMax,
                                                                                      Real_v &distance)
{
  // A fully planar GenericTrap is convex. On the host, clip directly against
  // its fixed lateral and z planes to avoid generic facet-container loads. The
  // compact facet representation remains faster on CUDA and is used below.
#ifndef VECCORE_CUDA
  if (unplaced.fTslHelper && (unplaced.fFlags & 0xF0u) == 0) {
    const Real_v invdirz = Real_v(1.) / NonZero(direction.z());
    const Real_v pz      = point.z();
    distance             = InfinityLength<Real_v>();
    if ((Abs(pz) - unplaced.fDz) > -kTolerance && pz * direction.z() >= Real_v(0.)) return;

    const Real_v ddz   = CopySign(unplaced.fDz, invdirz);
    const Real_v distz = -(pz + ddz) * invdirz;
    Real_v entry       = distz;
    Real_v limit       = Min(-(pz - ddz) * invdirz, stepMax);

    VECGEOM_PRAGMA_UNROLL(4)
    for (int i = 0; i < 4; ++i) {
      const Real_v projection = unplaced.fSurf[i].DotPlaneNormal(direction);
      const Real_v planeDist  = unplaced.fSurf[i].EvaluatePlanar(point);
      const bool outside      = planeDist > -kTolerance;
      if (outside) {
        if (projection > -kTolerance) return;
        entry = Max(entry, -planeDist / projection);
      } else if (projection > kTolerance) {
        limit = Min(limit, -planeDist / projection);
      }
    }
    if (entry >= limit - kTolerance) return;
    distance = entry;
    return;
  }
#endif

  // Fast path for planar faces via tessellated helper
  if (unplaced.fTslHelper) {
    Precision invdirz = 1. / NonZero(direction.z());
    // The upper flag bits identify collapsed faces, which are omitted from the
    // helper buffer. Otherwise a planar GenericTrap always has four facets.
    if ((unplaced.fFlags & 0xF0u) == 0)
      distance = unplaced.fTslHelper->DistanceToIn<false, 4>(point, direction, invdirz, stepMax);
    else
      distance = unplaced.fTslHelper->DistanceToIn<false>(point, direction, invdirz, stepMax);
    return;
  }
  // Summary:
  //  1) Handle trivially-inside start (convention distance=-1).
  //  2) Intersect the AABB to get an upper bound and early-out if unreachable.
  //  3) Try entering through Z slabs (cheap) and return if successful.
  //  4) Approact to bounding box distance w/o overshooting into. Improves numerical stability.
  //  5) Hoist v-window once, then test lateral faces with [lmin,lmax) window.
#ifndef VECGEOM_NO_WRONGSIDE_CONV
  bool completelyinside, completelyoutside;
  GenericKernelForContainsAndInside<Real_v, true>(unplaced, point, completelyinside, completelyoutside);

  // Already inside -> signal with -1 as convention
  if (completelyinside) {
    distance = Real_v(-1.);
    return;
  }
#endif

  // Intersect with bounding box first to get an upper bound
  Real_v bbdistance = Real_v(kInfLength);
  BoxImplementation::DistanceToIn(BoxStruct<Precision>(unplaced.fBBdimensions), point - unplaced.fBBorigin, direction,
                                  stepMax, bbdistance);

  // Initialize best distance as +inf; if we cannot reach the AABB we cannot hit the solid
  distance = InfinityLength<Real_v>();

  // If we cannot reach the AABB, we cannot reach the solid either
  bool done = bbdistance >= InfinityLength<Real_v>();
  if (done) return;

  // Potential hit with Z planes? Only if we are outside in Z and heading inward
  Real_v zsafety = Abs(point.z()) - unplaced.fDz;
  bool canhitz   = (zsafety > -kToleranceDist<Real_v>) && (point.z() * direction.z() < 0);

  if (canhitz) {
    // Distance to the relevant Z plane (as in Box algorithm)
    Real_v next = zsafety / Abs(direction.z());

    // Check containment on the Z plane polygon (top if dir.z>0, bottom otherwise)
    auto const vertices = (direction.z() > 0) ? &unplaced.fVertices[0] : &unplaced.fVertices[4];
    auto hitz           = InsideSection(point + next * direction, vertices);
    if (hitz) {
      distance = next;
      return;
    }
  }

  // ---------- Approach to (just outside) the AABB entry --------------
  // We start subsequent work from p0 = point + t_approach*dir, with t_approach = max(0, bbdistance - eps),
  Real_v t_approach   = Max(Real_v(0), bbdistance - kToleranceArb4<Real_v>);
  Vector3D<Real_v> p0 = point + t_approach * direction;
  Real_v step_budget  = stepMax - t_approach;
  if (step_budget <= Real_v(0)) {
    distance = InfinityLength<Real_v>();
    return;
  }

  // -------- Hoist v(t) in [0,1] window once (face-independent) --------
  // We gate by stepMax here to avoid generating candidates beyond the transport step.
  // IMPORTANT: We do not cap tmax by the AABB distance; instead we raise lmin later
  // with the (possibly finite) AABB entry so that faces cannot return hits before
  // the box is actually entered, while still allowing v-window to rule out ranges.
  Real_v tmin_v, tmax_v;
  if (!GenTrap_VInterval(p0.z(), direction.z(), unplaced.fDz2, step_budget, tmin_v, tmax_v)) {
    // Ray never reaches z within [-Dz,+Dz] in forward t within stepMax
    return;
  }

  // Do not accept hits before we actually pierce the AABB. From p0 the box entry is ~eps_bb ahead.
  Real_v lmin = Max(tmin_v, (t_approach > Real_v(0)) ? kToleranceArb4<Real_v> : Real_v(0));
  Real_v lmax = tmax_v;
  if (lmin > lmax) return;

  // Keep best local distance measured from p0; add t_approach at the end.
  Real_v best_local = InfinityLength<Real_v>();

  // Fall back to lateral surfaces and keep the minimum positive distance
  // Separate planar/twisted unrolled loops to reduce divergence on device compare to a mixed loop
  VECGEOM_PRAGMA_UNROLL(4)
  for (auto i = 0; i < 4; ++i) {
    if (!unplaced.IsDegenerated(i) && !unplaced.IsTwisted(i)) {
      // Planar faces: necessary front-facing check to avoid useless intersections
      Real_v dn = unplaced.fSurf[i].DotPlaneNormal(direction);
      if (dn >= -kToleranceDist<Real_v>) continue;
      // Cap face work by the z-feasible window upper bound
      // (we also pass lmin so the face can reject small distances early, before
      // any UV or normal work)
      Real_v face_limit = Min(best_local, lmax);
      best_local        = Min(best_local, DistanceToSurf(unplaced, p0, direction, i, false, lmin, face_limit));
    }
  }

  VECGEOM_PRAGMA_UNROLL(4)
  for (auto i = 0; i < 4; ++i) {
    if (!unplaced.IsDegenerated(i) && unplaced.IsTwisted(i)) {
      // Twisted faces: use the hoisted v-window to bound exact surface intersections
      // As above, we bound by tmax_v and provide lmin=tmin_v to allow
      // early rejection prior to UV/normal work.
      Real_v face_limit = Min(best_local, lmax);
      best_local        = Min(best_local, DistanceToSurf(unplaced, p0, direction, i, false, lmin, face_limit));
    }
  }

  // Report distance from the original point
  distance = (best_local >= InfinityLength<Real_v>()) ? InfinityLength<Real_v>() : (t_approach + best_local);
}

//______________________________________________________________________________
template <typename Real_v>
VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE void GenTrapImplementation::SetExitHit(UnplacedStruct_t const &unplaced,
                                                                                    ExitSurface surface,
                                                                                    SurfaceHitView<Real_v> *hitInfo)
{
  if (surface == ExitSurface::kNone) return;

  if (surface == ExitSurface::kZMin || surface == ExitSurface::kZMax) {
    hitInfo->SetLeftBehind(true);
    if (hitInfo->WantsNormal()) {
      const Real_v z = surface == ExitSurface::kZMin ? Real_v(-1.) : Real_v(1.);
      hitInfo->SetNormal(Vector3D<Real_v>(Real_v(0.), Real_v(0.), z));
    }
    return;
  }

  const int face = static_cast<int>(surface) - static_cast<int>(ExitSurface::kLateral0);
  // A planar lateral supports both convex end sections, hence every section
  // obtained by linear interpolation, even when adjacent laterals are twisted.
  const bool leftBehind = !unplaced.IsTwisted(face);
  hitInfo->SetLeftBehind(leftBehind);
  if (leftBehind && hitInfo->WantsNormal()) hitInfo->SetNormal(Vector3D<Real_v>(unplaced.fSurf[face].GetPlaneNormal()));
}

//______________________________________________________________________________
template <typename Real_v, bool RecordHitInfo>
VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE void GenTrapImplementation::DistanceToOutImpl(
    UnplacedStruct_t const &unplaced, Vector3D<Real_v> const &point, Vector3D<Real_v> const &direction,
    Real_v const & /*stepMax*/, Real_v &distance, SurfaceHitView<Real_v> *hitInfo)
{
  if constexpr (RecordHitInfo) hitInfo->Clear();

  if (unplaced.fTslHelper) {
    if constexpr (RecordHitInfo) {
      const Real_v safz = Abs(point.z()) - unplaced.fDz;
      if (safz > Real_v(kTolerance)) {
        distance = Real_v(-kTolerance);
        return;
      }

      const Real_v vz           = direction.z();
      const Real_v capDistance  = vz == Real_v(0.)
                                      ? InfinityLength<Real_v>()
                                      : Max((CopySign(Real_v(unplaced.fDz), vz) - point.z()) / vz, Real_v(0.));
      int face                  = -1;
      const Real_v sideDistance = unplaced.fTslHelper->DistanceToOut(point, direction, face);
      if (face < 0 && sideDistance < Real_v(0.)) {
        distance = sideDistance;
        return;
      }

      if (sideDistance < capDistance) {
        distance = sideDistance;
        hitInfo->SetLeftBehind(true);
        if (hitInfo->WantsNormal())
          hitInfo->SetNormal(Vector3D<Real_v>(unplaced.fTslHelper->GetFacet(face).GetNormal()));
      } else {
        distance = capDistance;
        SetExitHit(unplaced, vz < Real_v(0.) ? ExitSurface::kZMin : ExitSurface::kZMax, hitInfo);
      }
    } else {
      distance = unplaced.fTslHelper->DistanceToOut(point, direction);
    }
    return;
  }
  // With wrong-side conventions enabled, outside starts return -1. Otherwise
  // callers guarantee an inside/surface start. Compute
  // the minimum positive distance to any exiting surface (Z slabs + laterals).
  // Z slabs are attempted first as they are cheaper; lateral faces use the
  // same DistanceToSurf() with in=true and a shrinking lmax=distance.
  bool assumeInsideStart = true;
#ifndef VECGEOM_NO_WRONGSIDE_CONV
  distance = Real_v(-1.); // default if already outside
  bool completelyoutside;
  // Enforce the public wrong-side convention when requested. The classification
  // also identifies strict interior starts for the direct exit-root path.
  GenericKernelForContainsAndInside<Real_v, true>(unplaced, point, assumeInsideStart, completelyoutside);
  if (completelyoutside) return;
#endif

  // Candidate intersection with the Z planes (skip if nearly parallel)
  Real_v sign         = (direction.z() < 0) ? Real_v(-1.) : Real_v(1.);
  Real_v limit        = Real_v(1.); // <- This is a first implementation, we should use the real limit
  bool skipZ          = Abs(direction.z()) * limit < kToleranceDist<Real_v>;
  distance            = skipZ ? InfinityLength<Real_v>() : (sign * unplaced.fDz - point.z()) / NonZero(direction.z());
  ExitSurface surface = skipZ ? ExitSurface::kNone : (sign < Real_v(0.) ? ExitSurface::kZMin : ExitSurface::kZMax);

  // Fall back to lateral surfaces and keep the minimum positive distance
  // Separate planar/twisted unrolled loops to reduce divergence on device compare to a mixed loop
  VECGEOM_PRAGMA_UNROLL(4)
  for (auto i = 0; i < 4; ++i) {
    if (!unplaced.IsDegenerated(i) && !unplaced.IsTwisted(i)) {
      if constexpr (RecordHitInfo) {
        const Real_v candidate =
            DistanceToSurf(unplaced, point, direction, i, true, Real_v(0), distance, assumeInsideStart);
        if (candidate < distance) {
          distance = candidate;
          surface  = static_cast<ExitSurface>(static_cast<int>(ExitSurface::kLateral0) + i);
        }
      } else {
        distance =
            Min(distance, DistanceToSurf(unplaced, point, direction, i, true, Real_v(0), distance, assumeInsideStart));
      }
    }
  }

  VECGEOM_PRAGMA_UNROLL(4)
  for (auto i = 0; i < 4; ++i) {
    if (!unplaced.IsDegenerated(i) && unplaced.IsTwisted(i)) {
      if constexpr (RecordHitInfo) {
        const Real_v candidate =
            DistanceToSurf(unplaced, point, direction, i, true, Real_v(0), distance, assumeInsideStart);
        if (candidate < distance) {
          distance = candidate;
          surface  = static_cast<ExitSurface>(static_cast<int>(ExitSurface::kLateral0) + i);
        }
      } else {
        distance =
            Min(distance, DistanceToSurf(unplaced, point, direction, i, true, Real_v(0), distance, assumeInsideStart));
      }
    }
  }

  if constexpr (RecordHitInfo) SetExitHit(unplaced, surface, hitInfo);
}

//______________________________________________________________________________
template <typename Real_v>
VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE void GenTrapImplementation::DistanceToOut(
    UnplacedStruct_t const &unplaced, Vector3D<Real_v> const &point, Vector3D<Real_v> const &direction,
    Real_v const &stepMax, Real_v &distance, SurfaceHitView<Real_v> *hitInfo)
{
  if (hitInfo)
    DistanceToOutImpl<Real_v, true>(unplaced, point, direction, stepMax, distance, hitInfo);
  else
    DistanceToOutImpl<Real_v, false>(unplaced, point, direction, stepMax, distance, nullptr);
}

//______________________________________________________________________________
template <typename Real_v>
VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE void GenTrapImplementation::SafetyToIn(UnplacedStruct_t const &unplaced,
                                                                                    Vector3D<Real_v> const &point,
                                                                                    Real_v &safety)
{
  if (unplaced.fTslHelper) {
    safety = unplaced.fTslHelper->SafetyToIn(point);
    return;
  }
  // If outside the bounding box, use the box safety (cheap and conservative)
  bool inside;
  BoxImplementation::Contains(BoxStruct<Precision>(unplaced.fBBdimensions), point - unplaced.fBBorigin, inside);
  if (!inside) {
    // This is not optimal if top and bottom faces are not on top of each other
    BoxImplementation::SafetyToIn(BoxStruct<Precision>(unplaced.fBBdimensions), point - unplaced.fBBorigin, safety);
    return;
  }

  // Combine Z safety with lateral faces using SafetyArb4
  safety = Abs(point[2]) - unplaced.fDz;
  safety = SafetyArb4(unplaced, point, safety, false);
}

//______________________________________________________________________________
template <typename Real_v>
VECCORE_ATT_HOST_DEVICE void GenTrapImplementation::SafetyToOut(UnplacedStruct_t const &unplaced,
                                                                Vector3D<Real_v> const &point, Real_v &safety)
{
  // Dispatcher for SafetyToOut
  if (unplaced.fTslHelper) {
    safety = unplaced.fTslHelper->SafetyToOut(point);
    return;
  }
  // Generic implementation for SafetyToOut
  // Safety convention:
  //  - positive value is a *lower bound* to the distance to exit
  //  - for Z, we return signed distance to the nearest slab
  //  - if already outside along Z (safety >= -tol), return -1 as sentinel

  // Z safety is simple: distance to nearest top/bottom plane (signed)
  safety = Abs(point[2]) - unplaced.fDz;

  // Already outside on Z? return -1 to signal this convention
  if (safety > -kTolerance) {
    safety = Real_v(-1.);
    return;
  }

  // Aggregate with lateral faces
  safety = SafetyArb4(unplaced, point, safety, true);
}

//______________________________________________________________________________
template <typename Real_v>
VECCORE_ATT_HOST_DEVICE void GenTrapImplementation::NormalKernel(UnplacedStruct_t const &unplaced,
                                                                 Vector3D<Real_v> const &point,
                                                                 Vector3D<Real_v> &normal, bool &valid)
{
  // Computes the normal on a surface and returns it as a unit vector
  // Strategy:
  //   * Prefer Z planes when |z| is close to Dz.
  //   * Identify incident lateral faces on the section at z and use cached
  //     plane normals or implicit gradients, averaging multiple matches.
  //   * Fall back to the closest projected edge for unmatched points; twisted
  //     fallback faces require (u,v) in a tolerance-expanded [0,1]^2.
  valid = true;
  normal.Set(0.);

  // First, check proximity to the top/bottom planes
  Real_v safz = Abs(unplaced.fDz - Abs(point.z()));
  if (safz < kTolerance) {
    normal[2] = Sign(point.z());
    return;
  }

  // Resolve incident lateral faces from the section at z before falling back
  // to closest-edge projections. The implicit gradient avoids UV recovery.
  Vector2D<Real_v> vert[4];
  FillVerticesAtZ(unplaced, point.z(), vert);
  int surfaceCount    = 0;
  bool firstIsTwisted = false;
  if (Abs(point.z()) < unplaced.fDz) {
    for (int i = 0; i < 4; ++i) {
      if (unplaced.IsDegenerated(i)) continue;

      const int next         = (i + 1) % 4;
      const Real_v tolerance = unplaced.IsTwisted(i) ? kToleranceArb4<Real_v> : Real_v(kTolerance);
      const Real_v minx      = Min(vert[i].x(), vert[next].x()) - tolerance;
      const Real_v maxx      = Max(vert[i].x(), vert[next].x()) + tolerance;
      const Real_v miny      = Min(vert[i].y(), vert[next].y()) - tolerance;
      const Real_v maxy      = Max(vert[i].y(), vert[next].y()) + tolerance;
      if (!(point.x() > minx && point.x() < maxx && point.y() > miny && point.y() < maxy)) continue;

      const bool twisted = unplaced.IsTwisted(i);
      if (!twisted) {
        if (!(Abs(unplaced.fSurf[i].EvaluatePlanar(point)) < tolerance)) continue;
      } else {
        const Real_v dx    = vert[next].x() - vert[i].x();
        const Real_v dy    = vert[next].y() - vert[i].y();
        const Real_v dpx   = point.x() - vert[i].x();
        const Real_v dpy   = point.y() - vert[i].y();
        const Real_v side  = dx * dpy - dy * dpx;
        const Real_v lensq = dx * dx + dy * dy;
        if (!(side * side < tolerance * tolerance * lensq)) continue;
      }

      const auto faceNormal = twisted ? unplaced.fSurf[i].EvaluateGradient(point) : unplaced.fSurf[i].GetPlaneNormal();
      if (surfaceCount == 0) {
        normal         = faceNormal;
        firstIsTwisted = twisted;
      } else {
        // Only edge hits normalize inside the loop. A regular face hit needs
        // one normalization after selection, including across divergent lanes.
        if (surfaceCount == 1 && firstIsTwisted) normal.Normalize();
        normal += twisted ? faceNormal.Normalized() : faceNormal;
      }
      ++surfaceCount;
    }
  }
  if (surfaceCount != 0) {
    if (surfaceCount > 1 || firstIsTwisted) normal.Normalize();
    return;
  }

  int iseg;
  GetClosestEdge<Real_v>(point, vert, iseg);

  if (!unplaced.IsTwisted(iseg)) {
    // Planar face: use precomputed outward normal
    normal = unplaced.fSurf[iseg].GetPlaneNormal();
    return;
  }

  Real_v u, v;
  UNormal(unplaced, point, iseg, normal, u, v);
  normal.Normalize();

  // Consider points very near boundaries as valid (u,v within tolerance-extended [0,1])
  valid = (u > -kToleranceDist<Real_v>) && (u < 1. + kToleranceDist<Real_v>) && (v > -kToleranceDist<Real_v>) &&
          (v < 1. + kToleranceDist<Real_v>);
}

//______________________________________________________________________________
template <typename Real_v>
VECCORE_ATT_HOST_DEVICE void GenTrapImplementation::GetClosestEdge(Vector3D<Real_v> const &point,
                                                                   Vector2D<Real_v> const vertxy[4], int &iseg)
{
  // Get index of the edge of the quadrilateral represented by vert closest to point.
  iseg = 0;
  Real_v lsq, dx, dy, dpx, dpy, u;
  Real_v safe = InfinityLength<Real_v>();
  Real_v ssq  = InfinityLength<Real_v>();

  for (int i = 0; i < 4; ++i) {
    int j = (i + 1) % 4;
    dx    = vertxy[j].x() - vertxy[i].x();
    dy    = vertxy[j].y() - vertxy[i].y();
    dpx   = point.x() - vertxy[i].x();
    dpy   = point.y() - vertxy[i].y();
    lsq   = dx * dx + dy * dy;

    // Skip degenerated segments
    if (lsq < kToleranceDistSquared<Real_v>) continue;

    // Projection fraction of P onto segment [i,j]
    u = (dpx * dx + dpy * dy);
    if (u > lsq) {
      // Outside, closer to vertex j
      dpx = point.x() - vertxy[j].x();
      dpy = point.y() - vertxy[j].y();
    } else if (u >= 0) {
      // Projection inside the segment
      auto invlsq = Real_v(1.) / lsq;
      dpx -= u * dx * invlsq;
      dpy -= u * dy * invlsq;
    }

    ssq = dpx * dpx + dpy * dpy;
    if (ssq < safe) {
      safe = ssq;
      iseg = i;
    }
  }
}

//______________________________________________________________________________
template <typename Real_v>
VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE void GenTrapImplementation::UNormal(UnplacedStruct_t const &unplaced,
                                                                                 Vector3D<Real_v> const &point, int i,
                                                                                 Vector3D<Real_v> &unorm, Real_v &u,
                                                                                 Real_v &v)
{
  // Get the (u, v) parameters corresponding to the cartesian point
  auto uv = XYZtoUV(unplaced, point, i);
  u       = uv[0];
  v       = uv[1];

  // Tangent vectors of the bilinear patch
  // Use extended caches when available
#if GENTRAP_USE_HYBRID_COEFFS
  const auto &E0 = unplaced.fE0[i]; // C - A
  const auto &E1 = unplaced.fE1[i]; // D - B
  const auto &G0 = unplaced.fG0[i]; // B - A
  const auto &G1 = unplaced.fG1[i]; // D - C
  auto dPdu      = -(Real_v(1.) - v) * E0 - v * E1;
  auto dPdv      = (Real_v(1.) - u) * G0 + u * G1;
#else
  auto j = (i + 1) % 4;
  // Corner notation: face spans A->C on bottom, B->D on top
  auto const &A = unplaced.fVertices[i];
  auto const &C = unplaced.fVertices[j];
  auto const &B = unplaced.fVertices[i + 4];
  auto const &D = unplaced.fVertices[j + 4];
  auto dPdu     = -(Real_v(1.) - v) * (C - A) - v * (D - B);
  auto dPdv     = (Real_v(1.) - u) * (B - A) + u * (D - C);
#endif

  // Un-normalized outward normal
  unorm = dPdu.Cross(dPdv);
}

//______________________________________________________________________________
template <typename Real_v>
VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE bool GenTrapImplementation::InsideSection(Vector3D<Real_v> const &point,
                                                                                       Vector3D<Real_v> const v[4])
{
  bool inside = true;
  bool degen  = true;

  // Half-plane test around the quadrilateral; uses robust cross/length checks
  for (int i = 0; i < 4; ++i) {
    const auto vij = v[(i + 1) % 4] - v[i];
    auto len2      = vij.Mag2();
    if (len2 < kToleranceDistSquared<Real_v>) continue; // degenerate edge
    degen      = false;
    auto cross = (point - v[i]).CrossZ(vij);
    inside     = cross * cross >= kToleranceDistSquared<Real_v> * len2 ? (cross >= 0) : false;
    if (!inside) break; // early exit if outside any edge
  }
  if (degen) return false; // degenerate polygon is treated as outside
  return inside;
}

//______________________________________________________________________________
template <typename Real_v>
VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE void GenTrapImplementation::FillVerticesAtZ(
    UnplacedStruct_t const &unplaced, Real_v z, Vector2D<Real_v> vertices[4])
{
  // Interpolate XY vertices along generators to the specified z-plane
  for (int i = 0; i < 4; i++)
    vertices[i] = unplaced.fMidXY[i] + unplaced.fT[i] * z;
}

//______________________________________________________________________________
template <typename Real_v>
VECCORE_ATT_HOST_DEVICE Vector2D<Real_v> GenTrapImplementation::XYZtoUV(UnplacedStruct_t const &unplaced,
                                                                        Vector3D<Real_v> const &point, int i)
{
  // Parametrize v from z: bottom plane z=-Dz -> v=0, top plane z=+Dz -> v=1
  Real_v v = unplaced.fDz2 * (point[2] + unplaced.fDz);

  // Points E, F on the generators (A, C) and (B, D) at z = point[2]
#if GENTRAP_USE_HYBRID_COEFFS
  const auto &A  = unplaced.fVertices[i];
  const auto &G0 = unplaced.fG0[i];                // B - A
  const auto &E0 = unplaced.fE0[i];                // C - A
  const auto &E1 = unplaced.fE1[i];                // D - B
  auto E         = A + v * G0;                     // (1-v)A + vB
  auto EF        = (Real_v(1.) - v) * E0 + v * E1; // F - E but stable form
#else
  auto j  = (i + 1) % 4;
  auto E  = (1. - v) * unplaced.fVertices[i] + v * unplaced.fVertices[i + 4];
  auto F  = (1. - v) * unplaced.fVertices[j] + v * unplaced.fVertices[j + 4];
  auto EF = F - E;
#endif
  auto EP = point - E;

  // Project along EF to obtain u in [0,1]
  auto u = EP.Dot(EF) / NonZero(EF.Dot(EF));
  return Vector2D<Real_v>(u, v);
}

//______________________________________________________________________________
template <typename Real_v>
VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE bool GenTrapImplementation::InSurfLimits(UnplacedStruct_t const &unplaced,
                                                                                      Vector3D<Real_v> const &point,
                                                                                      int isurf)
{
  auto uv = XYZtoUV(unplaced, point, isurf);
  // Tolerant test against the unit square in (u,v). This admits points that
  // are numerically on the boundaries (within kToleranceDist) as "in".
  bool inside = (uv[0] > -kToleranceDist<Real_v>) && (uv[0] < Real_v(1.) + kToleranceDist<Real_v>) &&
                (uv[1] > -kToleranceDist<Real_v>) && (uv[1] < Real_v(1.) + kToleranceDist<Real_v>);
  return inside;
}

//______________________________________________________________________________
template <typename Real_v>
VECGEOM_FORCE_INLINE VECCORE_ATT_HOST_DEVICE Real_v GenTrapImplementation::SafetyArb4(UnplacedStruct_t const &unplaced,
                                                                                      Vector3D<Real_v> const &point,
                                                                                      Real_v safmax, bool inside)
{
  // See: https://gitlab.cern.ch/geant4/geant4-dev/-/merge_requests/5261
  Real_v safety = safmax;

  // Aggregate per-face conservative signed safeties
  VECGEOM_PRAGMA_UNROLL(4)
  for (int i = 0; i < 4; i++) {
    if (!unplaced.IsDegenerated(i) && !unplaced.IsTwisted(i)) {
      // Planar faces: signed distance to plane using stored outward normal
      Real_v sface = unplaced.fSurf[i].EvaluatePlanar(point);
      safety       = Max(safety, sface);
    }
  }

  VECGEOM_PRAGMA_UNROLL(4)
  // Twisted bilinear faces: use Lipschitz bound of the implicit function
  // Note: We could replace this for slightly twisted faces with the much cheaper
  //       but less precise: unplaced.fMesh[i].FastSignedSafety(point, inside);
  for (int i = 0; i < 4; i++) {
    if (!unplaced.IsDegenerated(i) && unplaced.IsTwisted(i)) {
      // Outside -> a face can only help if f(point) > 0
      if (!inside && unplaced.fSurf[i].Evaluate(point) <= Real_v(0)) continue;
      Real_v sface = unplaced.fSurf[i].SafetyLipschitz(point);
      safety       = Max(safety, sface);
    }
  }

  // Return signed conservative safety: negative when inside, positive outside.
  return (inside) ? -safety : safety;
}

} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif /* GENTRAPIMPLEMENTATION_H_ */
