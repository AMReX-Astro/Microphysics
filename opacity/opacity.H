#ifndef opacity_H
#define opacity_H

#include <actual_opacity.H>

AMREX_INLINE
void opacity_init ()
{
    actual_opacity_init();
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void opacity (Real& kp, amrex::Real& kr, amrex::Real rho,
              amrex::Real temp, amrex::Real rhoYe, amrex::Real nu,
              bool get_Planck_mean, bool get_Rosseland_mean)
{
    actual_opacity(kp, kr, rho, temp, rhoYe, nu, get_Planck_mean, get_Rosseland_mean);
}

#endif
