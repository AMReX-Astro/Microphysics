#ifndef MICROPHYSICS_MATH_H
#define MICROPHYSICS_MATH_H

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <esum.H>

// Override definitions of constants in GCEM to work around lack of support
// on Power9 for long double. We don't use arithmetic > 64 bit anyway.

#define GCEM_LOG_2 0.6931471805599453094172321214581765680755
#define GCEM_LOG_10 2.3025850929940456840179914546843642076011
#define GCEM_PI 3.1415926535897932384626433832795028841972
#define GCEM_LOG_PI 1.1447298858494001741434273513530587116473
#define GCEM_LOG_2PI 1.8378770664093454835606594728112352797228
#define GCEM_LOG_SQRT_2PI 0.9189385332046727417803297364056176398614
#define GCEM_SQRT_2 1.4142135623730950488016887242096980785697
#define GCEM_HALF_PI 1.5707963267948966192313216916397514420986
#define GCEM_SQRT_PI 1.7724538509055160272981674833411451827975
#define GCEM_SQRT_HALF_PI 1.2533141373155002512078826424055226265035
#define GCEM_E 2.7182818284590452353602874713526624977572

#include <gcem.hpp>

namespace microphysics_math {

    // quintic hermite polynomial functions
    // psi0 and its derivatives
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real psi0 (Real z)
    {
        return z * z * z * (z * (-6.0e0_rt * z + 15.0e0_rt) -10.0e0_rt) + 1.0e0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real dpsi0 (Real z)
    {
        return z * z * (z * (-30.0e0_rt*z + 60.0e0_rt) - 30.0e0_rt);
    }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real ddpsi0 (Real z)
    {
        return z * (z * (-120.0e0_rt * z + 180.0e0_rt) - 60.0e0_rt);
    }

    // psi1 and its derivatives
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real psi1 (Real z)
    {
        return z * ( z * z * (z * (-3.0e0_rt * z + 8.0e0_rt) - 6.0e0_rt) + 1.0e0_rt);
    }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real dpsi1 (Real z)
    {
        return z * z * ( z * (-15.0e0_rt * z + 32.0e0_rt) - 18.0e0_rt) + 1.0e0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real ddpsi1 (Real z)
    {
        return z * (z * (-60.0e0_rt * z + 96.0e0_rt) - 36.0e0_rt);
    }

    // psi2 and its derivatives
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real psi2 (Real z)
    {
        return 0.5e0_rt * z * z * ( z * ( z * (-z + 3.0e0_rt) - 3.0e0_rt) + 1.0e0_rt);
    }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real dpsi2 (Real z)
    {
        return 0.5e0_rt * z * (z * (z * (-5.0e0_rt * z + 12.0e0_rt) - 9.0e0_rt) + 2.0e0_rt);
    }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real ddpsi2 (Real z)
    {
        return 0.5e0_rt * (z * ( z * (-20.0e0_rt * z + 36.0e0_rt) - 18.0e0_rt) + 2.0e0_rt);
    }

    // cubic hermite polynomial functions
    // psi0 & derivatives
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real xpsi0(Real z)
    {
        return z * z * (2.0e0_rt * z - 3.0e0_rt) + 1.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real xdpsi0(Real z)
    {
        return z * (6.0e0_rt * z - 6.0e0_rt);
    }

    // psi1 & derivatives
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real xpsi1(Real z)
    {
        return z * (z * (z - 2.0e0_rt) + 1.0e0_rt);
    }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real xdpsi1(Real z)
    {
        return z * (3.0e0_rt * z - 4.0e0_rt) + 1.0e0_rt;
    }
}

#endif
