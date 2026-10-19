// Host stand-in for <loomcxx/scalar.h>, for compiling the f64 OCML recipes
// with a plain host compiler (tests/ocml/README.md).
//
// fmaf, truncf, absf, copysignf, sqrtf, floorf, roundevenf are exact IEEE
// operations, so the host results equal the Loom ones bit for bit.
// approximate::log2f/exp2f model v_log_f32/v_exp_f32 with the correctly
// rounded host functions; the hardware results can differ in the last bits.
// Only f64_math.H's cbrt uses them, and the pow/log10/tanh recipes do not.
#pragma once
#include <cmath>

namespace loom::scalar {
template <class F> inline F fmaf(F a, F b, F c) { return std::fma(a, b, c); }
template <class F> inline F truncf(F x) { return std::trunc(x); }
template <class F> inline F absf(F x) { return std::fabs(x); }
template <class F> inline F copysignf(F x, F y) { return std::copysign(x, y); }
template <class F> inline F sqrtf(F x) { return std::sqrt(x); }
template <class F> inline F floorf(F x) { return std::floor(x); }
template <class F> inline F roundevenf(F x) { return std::nearbyint(x); }
namespace approximate {
template <class F> inline F log2f(F x) { return std::log2(x); }
template <class F> inline F exp2f(F x) { return std::exp2(x); }
} // namespace approximate
} // namespace loom::scalar
