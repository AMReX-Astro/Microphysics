// Platform layer for the Loom burn kernel (see README.md).
//
// The same source compiles three ways:
//  - under loom-import-cxx (__loom__ defined): device math is bound to Loom
//    scalar operations or to the OCML transcriptions in f64_math.H and
//    f64_math_extra.H;
//  - with the host compiler (native validation): math calls exactly what the
//    stock Microphysics code calls (std:: functions), so results compare bit
//    for bit against the stock burner;
//  - with hipcc and -DMPL_HIP_STANDALONE (the standalone reproducer,
//    assets/standalone/): an ordinary HIP __global__ kernel. Math calls the
//    std:: functions, which hipcc lowers to its OCML like the stock HIP
//    build, or, with -DMPL_HIP_MATH_RECIPES, the same f64 recipes as the
//    Loom path (their <loomcxx/scalar.h> is then the device shim in
//    assets/f64-tests/ocml_gpu/device_shim/).
//    It is opt-in, not keyed on __HIPCC__: the test_react HIP build also
//    compiles these headers with hipcc (loom_launch.H) and needs the host
//    forms.
// MPL_GPU_KERNEL is defined for both GPU forms (Loom and standalone HIP).
#ifndef MPL_SUPPORT_H
#define MPL_SUPPORT_H

#if defined(MPL_HIP_MATH_RECIPES) && !defined(MPL_HIP_STANDALONE)
#error "MPL_HIP_MATH_RECIPES needs MPL_HIP_STANDALONE"
#endif

#if defined(__loom__)
#include <hip/hip_runtime.h>
#define MPL_GPU_KERNEL 1
#define MPL_DEVICE __device__
#define MPL_KERNEL __global__ [[loom::workgroup_size(128, 1, 1)]]
#define DEVICE __device__   // f64_math.H expects DEVICE
#include "f64_math.H"
#include "f64_math_extra.H"
#undef DEVICE
// A math function with no f64 recipe yet: referencing this symbol makes the
// Loom import fail with its name instead of silently using something else.
#define MPL_LOOM_TODO(name) mpl_loom_todo_missing_f64_recipe_##name
#elif defined(MPL_HIP_STANDALONE)
#include <hip/hip_runtime.h>
#include <cmath>
#define MPL_GPU_KERNEL 1
#define MPL_DEVICE __device__
#define MPL_KERNEL __global__ __launch_bounds__(128)
#if defined(MPL_HIP_MATH_RECIPES)
#define DEVICE __device__
#include "f64_math.H"
#include "f64_math_extra.H"
#undef DEVICE
#endif
#else
#include <cmath>
#define MPL_DEVICE inline
#define MPL_KERNEL inline
#endif

// Value-range facts for the Loom bounds checker (SUBRANGE/026 "access byte
// end may reach ..."): needed where a loop index into a private array has
// data-dependent bounds. No-op on the host, so validated host code is unchanged.
// Loom only accepts plain operands: compare LOCAL scalar variables with
// constants, e.g. MPL_ASSUME(1 <= j && j <= N - 1). A struct field
// (vs->L: "predicate member paths require direct non-static record fields")
// or arithmetic on the variable side (j + 1 <= N: "predicate operands must be
// retained scalar values or pure constants") fails import; copy the field into
// a local first and move constants to the right-hand side.
#if defined(__loom__)
#define MPL_ASSUME(cond) __builtin_assume(cond)
#else
#define MPL_ASSUME(cond) ((void)0)
#endif

namespace mpl {

using Real = double;

// amrex::Array1D / amrex::Array2D stand-ins with the same 1-based accessors,
// so translated code keeps the source's index expressions (a(i), m(i, j)).
// Loom imports them (probes o02, o11, o12, o15): pass them by reference, or
// reach them as members through a pointer (vs->y(n)); calling a member's
// accessor through a *struct reference* does not import, const or not
// (r19, r23).
// Storage is padded from index 0 (slot 0 unused for 1-based arrays) and
// Array2D is row-major, a(i, j) = arr[i * (YHI + 1) + j]; this is NOT AMReX's
// column-major layout, which no translated code depends on.
template <typename T, int XLO, int XHI>
struct Array1D {
    static_assert(0 <= XLO && XLO <= XHI);
    T arr[XHI + 1];
    MPL_DEVICE T& operator()(int i) { return arr[i]; }
    MPL_DEVICE const T& operator()(int i) const { return arr[i]; }
};

template <typename T, int XLO, int XHI, int YLO, int YHI>
struct Array2D {
    static_assert(0 <= XLO && XLO <= XHI && 0 <= YLO && YLO <= YHI);
    T arr[(XHI + 1) * (YHI + 1)];
    MPL_DEVICE T& operator()(int i, int j) { return arr[i * (YHI + 1) + j]; }
    MPL_DEVICE const T& operator()(int i, int j) const { return arr[i * (YHI + 1) + j]; }
};

namespace math {

// amrex::Math::powi (AMReX_Math.H): keep its multiplication grouping.
template <int N>
MPL_DEVICE Real powi(Real x) {
    if constexpr (N < 0) {
        return Real(1) / powi<-N>(x);
    } else if constexpr (N == 0) {
        return Real(1);
    } else if constexpr (N == 1) {
        return x;
    } else if constexpr (N == 2) {
        return x * x;
    } else if constexpr (N % 2 == 0) {
        return powi<2>(powi<N / 2>(x));
    } else {
        return x * powi<N - 1>(x);
    }
}

// amrex::min/max -> std::min/max: the FIRST operand wins on ties and on
// unordered comparisons. Not fmin/fmax (different NaN and signed-zero rules).
template <typename T>
MPL_DEVICE T min(T a, T b) { return (b < a) ? b : a; }
template <typename T>
MPL_DEVICE T max(T a, T b) { return (a < b) ? b : a; }

// amrex::Clamp (AMReX_Algorithm.H): (v < lo) ? lo : (hi < v) ? hi : v,
// returned by value instead of by const reference.
template <typename T>
MPL_DEVICE T clamp(T v, T lo, T hi) { return (v < lo) ? lo : ((hi < v) ? hi : v); }

#if defined(__loom__) || defined(MPL_HIP_MATH_RECIPES)
MPL_DEVICE Real abs(Real x) { return loom::scalar::absf(x); }
MPL_DEVICE Real fabs(Real x) { return loom::scalar::absf(x); }
MPL_DEVICE Real sqrt(Real x) { return loom::scalar::sqrtf(x); }
MPL_DEVICE Real floor(Real x) { return loom::scalar::floorf(x); }
MPL_DEVICE Real copysign(Real x, Real y) { return loom::scalar::copysignf(x, y); }
MPL_DEVICE Real exp(Real x) { return chemistry::f64_math::exp(x); }
MPL_DEVICE Real log(Real x) { return chemistry::f64_math::log(x); }
MPL_DEVICE Real cbrt(Real x) { return chemistry::f64_math::cbrt(x); }
MPL_DEVICE Real pow(Real x, Real y) { return chemistry::f64_math::pow(x, y); }
MPL_DEVICE Real log10(Real x) { return chemistry::f64_math::log10(x); }
MPL_DEVICE Real tanh(Real x) { return chemistry::f64_math::tanh(x); }
#else
MPL_DEVICE Real abs(Real x) { return std::abs(x); }
MPL_DEVICE Real fabs(Real x) { return std::fabs(x); }
MPL_DEVICE Real sqrt(Real x) { return std::sqrt(x); }
MPL_DEVICE Real floor(Real x) { return std::floor(x); }
MPL_DEVICE Real copysign(Real x, Real y) { return std::copysign(x, y); }
MPL_DEVICE Real exp(Real x) { return std::exp(x); }
MPL_DEVICE Real log(Real x) { return std::log(x); }
MPL_DEVICE Real cbrt(Real x) { return std::cbrt(x); }
MPL_DEVICE Real pow(Real x, Real y) { return std::pow(x, y); }
MPL_DEVICE Real log10(Real x) { return std::log10(x); }
MPL_DEVICE Real tanh(Real x) { return std::tanh(x); }
#endif

} // namespace math
} // namespace mpl

#endif
