#ifndef MICROPHYSICS_APPROX_MATH_H
#define MICROPHYSICS_APPROX_MATH_H

#include <AMReX_REAL.H>
#include <AMReX_GpuUtility.H>
#include <microphysics_math.H>
#include <cassert>

using namespace amrex::literals;

AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real fast_atan(const amrex::Real x) {
    ///
    /// Fast atan calculation for Positive x only!
    /// This is meant to be used in screening calculation,
    /// and the input, gamma, is always positive.
    /// Also note that the typical input is < 1.0 for light nuclei
    /// and ~ 10 for heavy nuclei at high density.
    /// So this function is optimized with this scenario.
    ///

    // Make sure input is positive

    assert(x > 0.0_rt);

    // define some constants

    constexpr amrex::Real PI_2 = GCEM_PI / 2.0_rt;
    constexpr amrex::Real PI_4 = GCEM_PI / 4.0_rt;

    // Return 0 if x = 0

    if (x == 0.0_rt) {
        return x;
    }

    // Check for nan

    if (amrex::isnan(x)) {
        return x;
    }

    // Check for infinity

    if (amrex::isinf(x) && x > 0.0_rt) {
        return PI_2;
    }

    // Check for ones

    if (x == 1.0_rt) {
        return PI_4;
    }

    // Now calculate Atan(x) using Pade Approx

    amrex::Real y{0.0_rt};

    if (x < 1.0_rt) {
        // formula 1;
        // y = PI_4*x - x*(x - 1.0_rt) * (0.2447_rt + 0.0663*x);

        // formula 2;

        constexpr amrex::Real A = 0.0776509570923569;
        constexpr amrex::Real B = -0.287434475393028;
        constexpr amrex::Real C = GCEM_PI / 4.0_rt - A - B;

        amrex::Real x2 = x*x;
        y = ((A*x2 + B)*x2 + C)*x;
    }

    // // positive case
    // if (x > 0.0_rt) {
    //     const amrex::Real x2 = x * x;
    //     const amrex::Real x3 = x2 * x;
    //     const amrex::Real x4 = x3 * x;
    //     const amrex::Real x6 = x4 * x2;

    //     amrex::Real num = 7.0_rt * (23.0_rt*(x6 - 1.0_rt) +
    //                                 33.0_rt*x*(x4 - 1.0_rt) +
    //                                 30.0_rt*x2*(x2 + 1.0_rt));

    //     amrex::Real denom = 5.0_rt * (41.0_rt*(x6 + 1.0_rt) +
    //                                   111.0_rt*x*(x4 + 1.0_rt) +
    //                                   195.0_rt*x2*(x2 + 1.0_rt) +
    //                                   230.0_rt*x3);

    //     y = num / denom + PI_4;
    // }

    return y;
}

#endif
