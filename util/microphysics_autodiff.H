#ifndef MICROPHYSICS_AUTODIFF_H
#define MICROPHYSICS_AUTODIFF_H

#include <type_traits>

#include <AMReX.H>
#include <AMReX_Algorithm.H>
#include <AMReX_REAL.H>
#include <AMReX_Array.H>

#include <approx_math.H>

// required for AMREX_GPU_HOST_DEVICE, which is used via AUTODIFF_DEVICE_FUNC
#include <AMReX_GpuQualifiers.H>
// disable some optimizations that break standard left-to-right operator
// associativity, giving slightly different results with Dual vs. double
#define AUTODIFF_STRICT_ASSOCIATIVITY
#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/utils/derivative.hpp>

namespace microphysics_autodiff {

/**
 * A static-sized array that supports the math operations needed on the gradient
 * component of a Dual number.
 *
 * \tparam XLO Index for lower bound. Can be other than 0.
 * \tparam XHI Index for upper bound.
 */
template<typename T, int XLO, int XHI>
struct GradArray : public amrex::Array1D<T, XLO, XHI>
{
    /**
     * Member types
     */
    using value_type = T;
    using size_type = unsigned int;


    /// Construct a GradArray whose components are zero.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    constexpr GradArray() noexcept : GradArray(0) {}

    /// Construct a GradArray whose components are all equal to a scalar.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    explicit constexpr GradArray(T s) noexcept : amrex::Array1D<T, XLO, XHI>() {
        for(size_type i = 0; i < this->size(); ++i) {
            this->arr[i] = s;
        }
    }

    /**
     * Arithmetic operators
     */

    /// Modifies this GradArray by component-wise addition by argument.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    constexpr GradArray& operator+=(const GradArray& p) noexcept
    {
        for(size_type i = 0; i < this->size(); ++i) {
            this->arr[i] += p.arr[i];
        }
        return *this;
    }

    /// Modifies this GradArray by component-wise subtraction by argument.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    constexpr GradArray& operator-=(const GradArray& p) noexcept
    {
        for(size_type i = 0; i < this->size(); ++i) {
            this->arr[i] -= p.arr[i];
        }
        return *this;
    }

    /// Modifies this GradArray by multiplying each component by a scalar.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    constexpr GradArray& operator*=(T s) noexcept
    {
        for(size_type i = 0; i < this->size(); ++i) {
            this->arr[i] *= s;
        }
        return *this;
    }

    /// Modifies this GradArray by component-wise multiplication by argument.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    constexpr GradArray& operator*=(const GradArray& p) noexcept
    {
        for(size_type i = 0; i < this->size(); ++i) {
            this->arr[i] *= p.arr[i];
        }
        return *this;
    }

    /// Modifies this GradArray by dividing each component by a scalar.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    constexpr GradArray& operator/=(T s) noexcept
    {
        for(size_type i = 0; i < this->size(); ++i) {
            this->arr[i] /= s;
        }
        return *this;
    }

    /// Modifies this GradArray by component-wise division by argument.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    constexpr GradArray& operator/=(const GradArray& p) noexcept
    {
        for(size_type i = 0; i < this->size(); ++i) {
            this->arr[i] /= p.arr[i];
        }
        return *this;
    }
};

/**
 * Unary operators
 */
template<typename T, int XLO, int XHI>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GradArray<T, XLO, XHI> operator-(const GradArray<T, XLO, XHI>& s) noexcept
{
    GradArray<T, XLO, XHI> tmp(s);
    tmp *= -1;
    return tmp;
}

/**
 * Binary arithmetic operators: GradArray @ scalar
 */
template<typename T, int XLO, int XHI>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GradArray<T, XLO, XHI> operator*(const GradArray<T, XLO, XHI>& p, T s) noexcept
{
    GradArray<T, XLO, XHI> tmp(p);
    tmp *= s;
    return tmp;
}

template<typename T, int XLO, int XHI>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GradArray<T, XLO, XHI> operator/(const GradArray<T, XLO, XHI>& p, T s) noexcept
{
    GradArray<T, XLO, XHI> tmp(p);
    tmp /= s;
    return tmp;
}

/**
 * Binary arithmetic operators: scalar @ GradArray
 */
template<typename T, int XLO, int XHI>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GradArray<T, XLO, XHI> operator*(T s, const GradArray<T, XLO, XHI>& p) noexcept
{
    GradArray<T, XLO, XHI> tmp(p);
    tmp *= s;
    return tmp;
}

template<typename T, int XLO, int XHI>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GradArray<T, XLO, XHI> operator/(T s, const GradArray<T, XLO, XHI>& p) noexcept
{
    GradArray<T, XLO, XHI> tmp(s);
    tmp /= p;
    return tmp;
}

/**
 * Binary arithmetic operators: GradArray @ GradArray
 */
template<typename T, int XLO, int XHI>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GradArray<T, XLO, XHI> operator+(const GradArray<T, XLO, XHI>& s,
                                 const GradArray<T, XLO, XHI>& p) noexcept
{
    GradArray<T, XLO, XHI> tmp(s);
    tmp += p;
    return tmp;
}

template<typename T, int XLO, int XHI>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GradArray<T, XLO, XHI> operator-(const GradArray<T, XLO, XHI>& s,
                                 const GradArray<T, XLO, XHI>& p) noexcept
{
    GradArray<T, XLO, XHI> tmp(s);
    tmp -= p;
    return tmp;
}

template<typename T, int XLO, int XHI>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GradArray<T, XLO, XHI> operator*(const GradArray<T, XLO, XHI>& s,
                                 const GradArray<T, XLO, XHI>& p) noexcept
{
    GradArray<T, XLO, XHI> tmp(s);
    tmp *= p;
    return tmp;
}

template<typename T, int XLO, int XHI>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GradArray<T, XLO, XHI> operator/(const GradArray<T, XLO, XHI>& s,
                                 const GradArray<T, XLO, XHI>& p) noexcept
{
    GradArray<T, XLO, XHI> tmp(s);
    tmp /= p;
    return tmp;
}

// get<i>(arr) is needed to be able to do structured binding declarations like:
//    auto [dfdx, dfdy, dfdz] = autodiff::derivative(result);
template<std::size_t I, typename T, int XLO, int XHI>
T get(const GradArray<T, XLO, XHI>& arr) {
    return arr.arr[I];
}

} // namespace microphysics_autodiff

// std::tuple_size and std::tuple_element specializations are also needed for
// structured binding
template<typename T, int XLO, int XHI>
struct std::tuple_size<microphysics_autodiff::GradArray<T, XLO, XHI>>
    : public std::integral_constant<std::size_t, microphysics_autodiff::GradArray<T, XLO, XHI>::size()> {};

template<std::size_t I, typename T, int XLO, int XHI>
struct std::tuple_element<I, microphysics_autodiff::GradArray<T, XLO, XHI>> { using type = T; };

// open the autodiff namespace so we can make our own changes
namespace autodiff {
namespace detail {

/// ArithmeticTraits for GradArray (the array should be treated like a single number)
template<typename T, int XLO, int XHI>
struct ArithmeticTraits<microphysics_autodiff::GradArray<T, XLO, XHI>>
{
    static constexpr bool isArithmetic = true;
};

// add a couple of missing math functions

// natural logarithm of 1+x (std::log1p)
using std::log1p;

struct Log1pOp {};

template<typename R>
using Log1pExpr = UnaryExpr<Log1pOp, R>;

template<typename R, Requires<isExpr<R>> = true>
AUTODIFF_DEVICE_FUNC constexpr auto log1p(R&& r) -> Log1pExpr<R> { return { std::forward<R>(r) }; }

template<typename T, typename G>
AUTODIFF_DEVICE_FUNC constexpr void apply(Dual<T, G>& self, Log1pOp)
{
    const T aux = One<T>() / (1.0 + self.val);
    self.val = log1p(self.val);
    self.grad *= aux;
}


// cube root (std::cbrt)
using std::cbrt;

struct CbrtOp {};

template<typename R>
using CbrtExpr = UnaryExpr<CbrtOp, R>;

template <typename R, Requires<isExpr<R>> = true>
AUTODIFF_DEVICE_FUNC constexpr auto cbrt(R&& r) -> CbrtExpr<R> { return { std::forward<R>(r) }; }

template<typename T, typename G>
AUTODIFF_DEVICE_FUNC constexpr void apply(Dual<T, G>& self, CbrtOp)
{
    self.val = cbrt(self.val);
    self.grad *= 1.0 / (3.0 * self.val * self.val);
}

// sincos from amrex::Math
// returns two values, so we can't define it as a UnaryExpr
template <typename R, Requires<isExpr<R>> = true, typename dual_t = DualType<R>>
AUTODIFF_DEVICE_FUNC AMREX_FORCE_INLINE
auto sincos(R&& expr) -> std::pair<dual_t, dual_t>
{
    // ensure any passed expressions are evaluated first
    dual_t x = std::forward<R>(expr);
    std::pair<dual_t, dual_t> result;
    auto [sin, cos] = amrex::Math::sincos(x.val);
    result.first.val = sin;
    result.first.grad = x.grad;
    result.first.grad *= cos;

    result.second.val = cos;
    result.second.grad = std::move(x).grad;
    result.second.grad *= -sin;
    return result;
}

// exponentiation by a compile-time integer (amrex::Math::powi)

// this needs to go in its own namespace, to avoid an ambiguous overload error
// due to ADL from the recursive call in amrex::Math::powi
namespace aux {
template <int Power, typename R>
AUTODIFF_DEVICE_FUNC AMREX_FORCE_INLINE
constexpr auto powi(R&& r) {
    return amrex::Math::powi<Power>(autodiff::eval(std::forward<R>(r)));
}
}


// custom functions from approx_math.H

// fast_atan
struct FastAtanOp {};

template<typename R>
using FastAtanExpr = UnaryExpr<FastAtanOp, R>;

template<typename R, Requires<isExpr<R>> = true>
AUTODIFF_DEVICE_FUNC constexpr auto fast_atan(R&& r) -> FastAtanExpr<R> { return { std::forward<R>(r) }; }

template<typename T, typename G>
AUTODIFF_DEVICE_FUNC constexpr void apply(Dual<T, G>& self, FastAtanOp)
{
    const T aux = One<T>() / (1.0 + self.val * self.val);
    self.val = ::fast_atan(self.val);
    self.grad *= aux;
}

// fast_exp
struct FastExpOp {};

template <typename R>
using FastExpExpr = UnaryExpr<FastExpOp, R>;

template<typename R, Requires<isExpr<R>> = true>
AUTODIFF_DEVICE_FUNC constexpr auto fast_exp(R&& r) -> FastExpExpr<R> { return { std::forward<R>(r) }; }

template<typename T, typename G>
AUTODIFF_DEVICE_FUNC constexpr void apply(Dual<T, G>& self, FastExpOp)
{
    self.val = ::fast_exp(self.val);
    self.grad *= self.val;
}

} // namespace detail

// Redefine dual to use amrex::Real instead of double
using dual = HigherOrderDual<1, amrex::Real>;
/// Dual number type that can calculate multiple derivatives in a single pass.
template <int XLO, int XHI>
using dual_array = Dual<amrex::Real, microphysics_autodiff::GradArray<amrex::Real, XLO, XHI>>;

/**
 * Helper function to seed each of the input variables when using dual_array.
 */
template <int XLO, int XHI, typename... Args,
          typename Enable = std::enable_if_t<(... && std::is_same_v<Args, dual_array<XLO, XHI>>)>>
AUTODIFF_DEVICE_FUNC void seed_array(dual_array<XLO, XHI>& first_arg, Args&... args) {
    // first_arg is needed to deduce XLO and XHI, while Enable checks that all
    // of Args... are dual_arrays of the right size
    static_assert(1 + sizeof...(args) == microphysics_autodiff::GradArray<amrex::Real, XLO, XHI>::size(),
                  "number of arguments to seed_array does not match number of derivatives");
    int i = XLO;
    first_arg.grad(i) = 1.0_rt;
    ++i;
    // use a fold expression and immediately-invoked lambda to iterate over the
    // rest of the variadic arguments
    ([&] {
        args.grad(i) = 1.0_rt;
        ++i;
    } (), ...);
}

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
using amrex::Math::sincos, autodiff::detail::sincos;
using autodiff::detail::aux::powi;

using amrex::min, autodiff::detail::min;
using amrex::max, autodiff::detail::max;

using ::fast_atan, autodiff::detail::fast_atan;
using ::fast_exp, autodiff::detail::fast_exp;

} // namespace math_functions

} // namespace autodiff

namespace admath = autodiff::math_functions;

#endif
