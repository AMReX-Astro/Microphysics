#ifndef MICROPHYSICS_APPROX_MATH_H
#define MICROPHYSICS_APPROX_MATH_H

#include <AMReX_REAL.H>
#include <AMReX_GpuUtility.H>
#include <microphysics_math.H>
#include <cassert>
#include <cmath>

using namespace amrex::literals;

AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real atanf_1(const amrex::Real x) {
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
amrex::Real atanf(const amrex::Real x) {
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

    // Check for large number, close to infinity.
    // Error is ~1e-8 rad by not checking actual inf

    if (x > 1.8e8_rt) {
        return PI_2;
    }
    else if (x < -1.e8_rt) {
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
        y = PI_2 - atanf_1(1.0_rt / x);
    }
    else if (x < -1.0_rt) {
        y = -PI_2 - atanf_1(1.0_rt / x);
    }
    else {
        y = atanf_1(x);
    }

    return y;
}

#endif
