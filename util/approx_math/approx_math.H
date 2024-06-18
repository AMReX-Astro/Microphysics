#ifndef MICROPHYSICS_APPROX_MATH_H
#define MICROPHYSICS_APPROX_MATH_H

#include <AMReX_REAL.H>
#include <AMReX_GpuUtility.H>
#include <microphysics_math.H>
#include <cassert>
#include <cmath>

using namespace amrex::literals;

AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real fast_atan_1(const amrex::Real x) {
    ///
    /// This calculates atan within [-1, 1] range.
    ///

    // formula 1
    // constexpr amrex::Real PI_4 = GCEM_PI / 4.0_rt;
    // return PI_4 * x;

    // formula 2
    // constexpr amrex::Real A = 0.0776509570923569_rt;
    // constexpr amrex::Real B = -0.287434475393028_rt;
    // constexpr amrex::Real C = GCEM_PI / 4.0_rt - A - B;

    // amrex::Real x2 = x*x;
    // return ((A*x2 + B)*x2 + C)*x;

    // formula 3
    constexpr amrex::Real PI_4 = GCEM_PI / 4.0_rt;

    return PI_4*x - x*(std::abs(x) - 1.0_rt) *
        (0.2447_rt + 0.0663_rt*std::abs(x));

    // formula 4
    // constexpr amrex::Real PI_4 = GCEM_PI / 4.0_rt;
    // const amrex::Real x2 = x * x;
    // const amrex::Real x3 = x2 * x;
    // const amrex::Real x4 = x3 * x;
    // const amrex::Real x6 = x4 * x2;

    // amrex::Real num = 7.0_rt * (23.0_rt*(x6 - 1.0_rt) +
    //                             33.0_rt*x*(x4 - 1.0_rt) +
    //                             30.0_rt*x2*(x2 + 1.0_rt));

    // amrex::Real denom = 5.0_rt * (41.0_rt*(x6 + 1.0_rt) +
    //                               111.0_rt*x*(x4 + 1.0_rt) +
    //                               195.0_rt*x2*(x2 + 1.0_rt) +
    //                               230.0_rt*x3);

    // return num / denom + PI_4;

}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real fast_atan(const amrex::Real x) {
    ///
    /// Fast atan approximation calculations.
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
    else if (amrex::isinf(x) && x < 0.0_rt) {
        return -PI_2;
    }

    // Check for ones

    if (x == 1.0_rt) {
        return PI_4;
    }
    else if (x == -1.0_rt) {
        return -PI_4;
    }

    // Now calculate Atan(x) using approximations

    amrex::Real y;

    if (x > 1.0_rt) {
        y = PI_2 - fast_atan_1(1.0_rt / x);
    }
    else if (x < -1.0_rt) {
        y = -PI_2 - fast_atan_1(1.0_rt / x);
    }
    else {
        y = fast_atan_1(x);
    }

    return y;
}

#endif
