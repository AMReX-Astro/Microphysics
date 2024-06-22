#ifndef APPROX_MATH_H
#define APPROX_MATH_H

#include <AMReX_REAL.H>
#include <microphysics_math.H>

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

#endif
