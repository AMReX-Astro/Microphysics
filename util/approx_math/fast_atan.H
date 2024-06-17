#ifndef MICROPHYSICS_APPROX_MATH_H
#define MICROPHYSICS_APPROX_MATH_H

#include <AMReX_REAL.H>
#include <AMReX_GpuUtility.H>
#include <microphysics_math.H>

using namespace amrex::literals;

AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real fast_atan(const amrex::Real x) {
    ///
    /// This uses Pade-Approximation expanding around 1 and -1
    /// with order (m, n) = (6, 6).
    /// It gives incredible accuracy at all ranges.
    ///

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
    if (amrex::isinf(x) && x < 0.0_rt) {
        return -PI_2;
    }

    // Check for ones

    if (x == 1.0_rt) {
        return PI_4;
    }
    if (x == -1.0_rt) {
        return -PI_4;
    }

    // Now calculate Atan(x) using Pade Approx

    const amrex::Real x2 = x * x;
    const amrex::Real x3 = x2 * x;
    const amrex::Real x4 = x3 * x;
    const amrex::Real x5 = x4 * x;
    // const amrex::Real x6 = x5 * x;

    amrex::Real y{0.0_rt};

    // positive case
    if (x > 0.0_rt) {
        amrex::Real num = 7.0_rt * (23.0_rt*x6 + 33.0_rt*x5 + 30.0_rt*x4 -
                                    30.0_rt*x2 - 33.0_rt*x - 23.0_rt);
        amrex::Real denom = 5.0_rt * (41.0_rt*x6 + 111.0_rt*x5 + 195.0_rt*x4 +
                                      230.0_rt*x3 + 195.0_rt*x2 + 111.0_rt*x +
                                      41.0_rt);

        // amrex::Real num = 436.0_rt*x5 + 445.0_rt*x4 + 265.0_rt*x3 -
        //                   265.0_rt*x2 - 445.0_rt*x - 436.0_rt;
        // amrex::Real denom = 15.0_rt * (37.0_rt*x5 + 85.0_rt*x4 + 130.0_rt*x3 +
        //                                130.0_rt*x2 + 85.0_rt*x + 37.0_rt);
        y = num / denom + PI_4;
    }

    // negative case

    else {
        amrex::Real num = 7.0_rt * (-23.0_rt*x6 + 33.0_rt*x5 - 30.0_rt*x4 +
                                    30.0_rt*x2 - 33.0_rt*x + 23.0_rt);
        amrex::Real denom = 5.0_rt * (41.0_rt*x6 - 111.0_rt*x5 + 195.0_rt*x4 -
                                      230.0_rt*x3 + 195.0_rt*x2 - 111.0_rt*x +
                                      41.0_rt);
        // amrex::Real num = -436.0_rt*x5 + 445.0_rt*x4 - 265.0_rt*x3 -
        //                    265.0_rt*x2 + 445.0_rt*x - 436.0_rt;
        // amrex::Real denom = 15.0_rt * (37.0_rt*x5 - 85.0_rt*x4 + 130.0_rt*x3 -
        //                                130.0_rt*x2 + 85.0_rt*x - 37.0_rt);
        y = num / denom - PI_4;
    }

    return y;
}

#endif
