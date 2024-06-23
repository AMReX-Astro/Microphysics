#ifndef APPROX_MATH_H
#define APPROX_MATH_H

#include <AMReX_REAL.H>
#include <microphysics_math.H>
#include <cstdint>
#include <cstring>
#include <limits>

using namespace amrex::literals;

AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real fast_atan_1(const amrex::Real x) {
    ///
    /// This calculates atan within [-1, 1] range.
    ///

    ///
    /// Approximation of atan for x in [-1, 1]
    /// Max absolute error for this is ~0.0015 rad
    /// See Ref:
    /// Efficient Approximations for the Arctangent Function by Rajan 2006
    ///

    // constexpr amrex::Real PI_4 = GCEM_PI / 4.0_rt;

    // return PI_4*x - x*(std::abs(x) - 1.0_rt) *
    //     (0.2447_rt + 0.0663_rt*std::abs(x));


    ///
    /// Another approximation of atan for x in [-1, 1]
    /// Max absolute error for this is ~0.00063 rad
    /// See Ref:
    /// https://stackoverflow.com/questions/42537957/fast-accurate-atan-arctan-approximation-algorithm
    ///

    constexpr amrex::Real A = 0.0776509570923569_rt;
    constexpr amrex::Real B = -0.287434475393028_rt;
    constexpr amrex::Real C = GCEM_PI / 4.0_rt - A - B;

    amrex::Real x2 = x*x;
    return ((A*x2 + B)*x2 + C)*x;
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real fast_atan(const amrex::Real x) {
    ///
    /// Fast atan approximation calculations.
    ///

    constexpr amrex::Real PI_2 = 0.5_rt * GCEM_PI;

    ///
    /// If x < 0.113, then using arctan(x) ~ x
    /// gives better answer than the approximation below.
    /// And accuracy increase as x << 0.113.
    /// Also significantly faster.
    ///

    if (std::abs(x) < 0.113_rt) {
        return x;
    }

    // Check for large number, close to infinity.
    // Error is ~1e-8 rad by not checking actual inf

    if (x > 1.e8_rt) {
        return PI_2;
    }
    if (x < -1.e8_rt) {
        return -PI_2;
    }

    // Now calculate Atan(x) using approximations

    if (x > 1.0_rt) {
        return PI_2 - fast_atan_1(1.0_rt / x);
    }
    if (x < -1.0_rt) {
        return -PI_2 - fast_atan_1(1.0_rt / x);
    }
    return fast_atan_1(x);
}


///
/// A fast implementation of exp with double precision input
/// This gives reasonable accuracy across all range, ~2% error
///
/// This version uses memcpy to avoid potential undefined behavior
/// from type punning through a union in Ref (1) and (2)
///
/// Code is obtained from Ref (3):
/// Ref:
/// 1) A Fast, Compact Approximation of the Exponential Function
///    by Schraudolph 1999
/// 2) On a Fast Compact Approximation of the Exponential Function
///    by Cawley 2000
/// 3) https://gist.github.com/jrade/293a73f89dfef51da6522428c857802d
///

AMREX_GPU_HOST_DEVICE AMREX_INLINE
float jrade_exp(const float x) {
    /// For single precision input.
    /// a = 2^23 / ln2
    ///
    /// b = 2^23 * (x0 - C), where x0 = 127 is the exponent bias
    /// C = (ln(ln2 + 2/e) - ln2 - ln(ln2)) / ln2
    /// This is a constant shift term chosen to minimize maximum relative error
    ///
    /// Let C = ln(3/(8*ln2) + 0.5)/ln2
    /// in order to minimize RMS relative error.
    ///

    constexpr float a = gcem::pow(2.0F, 23) / 0.6931471805599453F;
    constexpr float b = gcem::pow(2.0F, 23) * (127.0F - 0.043677448F);
    float y = a * x + b;

    //
    // Return 0 for large negative number
    // Return Inf for large positive number
    //

    constexpr float c = gcem::pow(2.0F, 23);
    constexpr float d = gcem::pow(2.0F, 23) * 255.0F;
    if (y < c || y > d) {
        y = (y < c) ? 0.0F : d;
    }

    auto n = static_cast<std::int32_t>(y);
    memcpy(&y, &n, 4);
    return y;
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
double jrade_exp(const double x) {
    /// For double precision input.
    /// a = 2^23 / ln2
    ///
    /// b = 2^23 * (x0 - C), where x0 = 127 is the exponent bias
    /// C = (ln(ln2 + 2/e) - ln2 - ln(ln2)) / ln2
    /// This is a constant shift term chosen to minimize maximum relative error
    ///
    /// Let C = ln(3/(8*ln2) + 0.5)/ln2
    /// in order to minimize RMS relative error.
    ///

    constexpr double a = gcem::pow(2.0, 52) / 0.6931471805599453;
    constexpr double b = gcem::pow(2.0, 52) * (1023.0 - 0.04367744890362246);
    double y = a * x + b;

    //
    // Return 0 for large negative number
    // Return Inf for large positive number
    //

    constexpr double c = gcem::pow(2.0, 52);
    constexpr double d = gcem::pow(2.0, 52) * 2047.0;
    if (y < c || y > d) {
        y = (y < c) ? 0.0 : d;
    }

    auto n = static_cast<std::int64_t>(y);
    memcpy(&y, &n, 8);
    return y;
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
double ekmett_exp(const double x) {
    ///
    /// This is a more accurate than jrade_exp
    /// but it is roughly twice as slow.
    ///
    /// This uses the identity exp(x) = exp(x/2) / exp(-x/2),
    /// so there is extra factor of 0.5 in a
    ///
    /// This comes from:
    /// https://github.com/ekmett/approximate/blob/master/cbits/fast.c
    ///

    constexpr double a = gcem::pow(2.0, 52) / 0.6931471805599453;

    // For minimizes max relative error
    constexpr double b = gcem::pow(2.0, 52) * (1023.0 - 0.04367744890362246);

    // For minimizes rms relative error
    // constexpr double b = 0x3fdf127e83d16f12;

    double u = a * x + b;
    double v = b - a * x;

    //
    // Return 0 for large negative number
    // Return Inf for large positive number
    //

    constexpr double c = gcem::pow(2.0, 52);
    constexpr double d = gcem::pow(2.0, 52) * 2047.0;
    if (u < c || u > d) {
        u = (u < c) ? 0.0 : d;
    }

    auto n = static_cast<std::int64_t>(u);
    auto m = static_cast<std::int64_t>(v);

    memcpy(&u, &n, 8);
    memcpy(&v, &m, 8);

    return u / v;
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real fast_exp(const amrex::Real x) {
    ///
    /// Implementation of fast exponential.
    /// This combines Taylor series when x < 0.1
    /// and the fast exponential function algorithm from various sources:
    /// jrade: https://gist.github.com/jrade/293a73f89dfef51da6522428c857802d
    /// ekmett: https://github.com/ekmett/approximate/blob/master/cbits/fast.c
    ///

    // Use Taylor if number is smaller than 0.1
    // Minor performance hit, but much better accuracy when x < 0.1

    if (std::abs(x) < 0.1_rt) {
        return 1.0_rt + x;
    }

    return jrade_exp(x);
    // return ekmett_exp(x);
}
#endif
