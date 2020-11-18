#ifndef actual_opacity_H
#define actual_opacity_H

#include <AMReX.H>

#include <extern_parameters.H>

AMREX_INLINE
void actual_opacity_init () {}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_opacity (Real& kp, Real& kr, Real rho, Real temp, Real rhoYe, Real nu,
                     bool get_Planck_mean, bool get_Rosseland_mean)
{
    const Real Ksc = 0.4e0; // Thomson scattering
    const Real fac = 1.e-4; // Planck mean is assumed to be fac * Ksc

    if (get_Planck_mean) {
        kp = rhoYe * Ksc * fac;
    }

    if (get_Rosseland_mean) {
        kr = rhoYe * Ksc;
    }
}

#endif


