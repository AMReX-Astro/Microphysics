#ifndef NET_RATE_TYPE_H
#define NET_RATE_TYPE_H

#include <actual_network.H>

struct rate_t {
    amrex::Array1D<amrex::Real, 1, Rates::NumRates>  screened_rates;
#ifdef SCREENING
    amrex::Array1D<amrex::Real, 1, Rates::NumScreenPairs>  log_screen;
#endif
    amrex::Real enuc_weak;
};

struct rate_derivs_t {
    amrex::Array1D<amrex::Real, 1, Rates::NumRates>  screened_rates;
    amrex::Array1D<amrex::Real, 1, Rates::NumRates>  dscreened_rates_dT;
#ifdef SCREENING
    amrex::Array1D<amrex::Real, 1, Rates::NumScreenPairs>  log_screen;
    amrex::Array1D<amrex::Real, 1, Rates::NumScreenPairs>  dlog_screen_dT;
#endif
    amrex::Real enuc_weak;
    amrex::Array1D<amrex::Real, 1, NumSpec> denuc_weak_dY;
    amrex::Real denuc_weak_dT;
    amrex::Real drate_Fe52_n_n_to_Fe54_approx_dYN{};
    amrex::Real drate_Fe54_to_Fe52_n_n_approx_dYN{};
    amrex::Real drate_Fe54_n_n_to_Fe56_approx_dYN{};
    amrex::Real drate_Fe56_to_Fe54_n_n_approx_dYN{};
    amrex::Real drate_Ni56_n_n_to_Ni58_approx_dYN{};
    amrex::Real drate_Ni58_to_Ni56_n_n_approx_dYN{};
};


#endif