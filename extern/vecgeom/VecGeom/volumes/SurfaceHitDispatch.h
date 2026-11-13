/// @file SurfaceHitDispatch.h
/// @brief Compile-time detection for optional surface-hit-aware kernels.

#ifndef VECGEOM_VOLUMES_SURFACEHITDISPATCH_H_
#define VECGEOM_VOLUMES_SURFACEHITDISPATCH_H_

#include "VecGeom/base/Global.h"
#include "VecGeom/base/Vector3D.h"
#include "VecGeom/volumes/SurfaceHitView.h"

#include <type_traits>
#include <utility>

namespace vecgeom {

inline namespace VECGEOM_IMPL_NAMESPACE {

namespace SurfaceHitDispatch {

template <typename Impl, typename Struct, typename Real_t, typename = void>
struct HasDistanceToOutWithHit : std::false_type {};

template <typename Impl, typename Struct, typename Real_t>
struct HasDistanceToOutWithHit<
    Impl, Struct, Real_t,
    std::void_t<decltype(Impl::DistanceToOut(std::declval<Struct const &>(), std::declval<Vector3D<Real_t> const &>(),
                                             std::declval<Vector3D<Real_t> const &>(), std::declval<Real_t const &>(),
                                             std::declval<Real_t &>(), std::declval<SurfaceHitView<Real_t> *>()))>>
    : std::true_type {};

} // namespace SurfaceHitDispatch

} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif // VECGEOM_VOLUMES_SURFACEHITDISPATCH_H_
