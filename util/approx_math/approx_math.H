#ifndef APPROX_MATH_H
#define APPROX_MATH_H

#include <AMReX_REAL.H>
#include <microphysics_math.H>
#include <cstdint>
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


AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real cawley_exp(const double x) {
    ///
    /// A fast implementation of exp.
    /// Roughly 6-times faster than std::exp.
    /// With reasonable accuracy across all range.
    ///
    /// Code is obtained from Ref (2):
    /// Ref:
    /// 1) A Fast, Compact Approximation of the Exponential Function
    ///    by Schraudolph 1999
    /// 2) On a Fast Compact Approximation of the Exponential Function
    ///    by Cawley 2000
    ///


    //
    // If we use two 4 byte integer, then A = 2^20 / ln2
    // Because it is now represented by a single 8 byte integer,
    // we need to shift it by 4 additional bytes so 20+4*8=52
    //

    constexpr std::int64_t A = 6497320848556798;

    //
    // For two 4-byte integer representation, this is
    // 2^20*1023 - C, where 1023 is the exponent bias x0 and
    // C=60801 is the constant shift that minimizes error.
    //
    // For a single 8-byte integer, this is then
    // 2^52*1023 - C, where C is again the shift which is
    // calculated as 2^52*ln(3/(8*ln2) + 0.5)/ln2
    // The overall result is represented by B
    // And it is also an integer unlike two 4-byte integer case.
    //

    constexpr std::int64_t B = 0x3fef127e83d16f12LL;

    union {
        double d;
        std::int64_t i;
    } _eco;

    _eco.i = static_cast<std::int64_t>(A * x) + B;
    return _eco.d;
}



AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real fast_exp(const amrex::Real x) {
    ///
    /// Implementation of fast exponential.
    /// This combines Taylor series when x < 0.1
    /// and the fast exponential function algorithm
    /// described in Cawley 2000
    ///

    // Use Taylor if number is smaller than 0.1
    // It is ~1/10th slower than simply calling _fast_exp
    // if this is not being used due to the if-statement overhead.
    // But _fast_exp gives poor accuracy for x < 0.1

    if (std::abs(x) < 0.1_rt) {
        return 1.0_rt + x;
    }

    //
    // Cawley_exp won't handle large positive
    // or negative x properly.
    //

    if (x > 700.0_rt) {
        return std::numeric_limits<amrex::Real>::infinity();
    }

    if (x < -700.0_rt) {
        return 0.0_rt;
    }

    return cawley_exp(x);
}
#endif
