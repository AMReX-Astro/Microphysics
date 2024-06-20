#ifndef MICROPHYSICS_AUTODIFF_H
#define MICROPHYSICS_AUTODIFF_H

#include <AMReX.H>
#include <AMReX_Algorithm.H>
#include <AMReX_REAL.H>

// required for AMREX_GPU_HOST_DEVICE, which is used via AUTODIFF_DEVICE_FUNC
#include <AMReX_GpuQualifiers.H>
#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/utils/derivative.hpp>

// open the autodiff namespace so we can make our own changes
namespace autodiff {
namespace detail {

// add a couple of missing math functions

// natural logarithm of 1+x (std::log1p)
struct Log1pOp {};

template<typename R>
using Log1pExpr = UnaryExpr<Log1pOp, R>;

template<typename R, Requires<isExpr<R>> = true>
constexpr auto log1p(R&& r) -> Log1pExpr<R> { return { std::forward<R>(r) }; }

template<typename T, typename G>
constexpr void apply(Dual<T, G>& self, Log1pOp)
{
    const T aux = One<T>() / (1.0 + self.val);
    self.val = std::log1p(self.val);
    self.grad *= aux;
}


// cube root (std::cbrt)
struct CbrtOp {};

template<typename R>
using CbrtExpr = UnaryExpr<CbrtOp, R>;

template <typename R, Requires<isExpr<R>> = true>
constexpr auto cbrt(R&& r) -> CbrtExpr<R> { return { std::forward<R>(r) }; }

template<typename T, typename G>
constexpr void apply(Dual<T, G>& self, CbrtOp)
{
    self.val = std::cbrt(self.val);
    self.grad *= 1.0 / (3.0 * self.val * self.val);
}

} // namespace detail

// Redefine dual to use amrex::Real instead of double
using dual = HigherOrderDual<1, amrex::Real>;

// A new namespace that has both the STL math functions and the overloads for
// dual numbers, so we can write the same function name whether we're operating
// on autodiff::dual or amrex::Real.
namespace math_functions {

using std::abs, autodiff::detail::abs;
using std::acos, autodiff::detail::acos;
using std::asin, autodiff::detail::asin;
using std::atan, autodiff::detail::atan;
using std::atan2, autodiff::detail::atan2;
using std::cos, autodiff::detail::cos;
using std::exp, autodiff::detail::exp;
using std::log10, autodiff::detail::log10;
using std::log, autodiff::detail::log;
using std::pow, autodiff::detail::pow;
using std::sin, autodiff::detail::sin;
using std::sqrt, autodiff::detail::sqrt;
using std::tan, autodiff::detail::tan;
using std::cosh, autodiff::detail::cosh;
using std::sinh, autodiff::detail::sinh;
using std::tanh, autodiff::detail::tanh;
using std::erf, autodiff::detail::erf;
using std::hypot, autodiff::detail::hypot;

using std::log1p, autodiff::detail::log1p;
using std::cbrt, autodiff::detail::cbrt;

using amrex::min, autodiff::detail::min;
using amrex::max, autodiff::detail::max;

} // namespace math_functions

} // namespace autodiff

namespace admath = autodiff::math_functions;

#endif
