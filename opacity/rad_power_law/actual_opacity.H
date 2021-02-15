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
    Real nup_kpp = std::pow(nu, kappa_p_exp_p);
    Real nup_kpr = std::pow(nu, kappa_r_exp_p);
    Real nup_kps = std::pow(nu, scatter_exp_p);

    Real teff = amrex::max(temp, 1.0e-50_rt);
    teff = teff + rad_temp_floor * std::exp(-teff / (rad_temp_floor + 1.e-50_rt));

    Real ks = const_scatter * std::pow(rho, scatter_exp_m) * std::pow(teff, -scatter_exp_n) * nup_kps;

    if (get_Planck_mean) {
#ifndef AMREX_USE_GPU
        if (const_kappa_p < 0.0_rt) {
            amrex::Error("Must set Planck opacity constant");
        }
#endif
        kp = const_kappa_p * std::pow(rho, kappa_p_exp_m) * std::pow(teff, -kappa_p_exp_n) * nup_kpp;
        kp = amrex::max(kp, kappa_floor);
    }

    if (get_Rosseland_mean) {
#ifndef AMREX_USE_GPU
        if (const_kappa_r < 0.0_rt) {
            amrex::Error("Must set Rosseland opacity constant");
        }
#endif
        kr = const_kappa_r * std::pow(rho, kappa_r_exp_m) * std::pow(teff, -kappa_r_exp_n) * nup_kpr;
        kr = amrex::max(kr + ks, kappa_floor);
    }
}

#endif
