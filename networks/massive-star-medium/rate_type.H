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
    amrex::Array1D<amrex::Real, 1, NumSpec>  dweak_ydot_dYe;
    amrex::Real enuc_weak;
    amrex::Array1D<amrex::Real, 1, NumSpec> denuc_weak_dY;
    amrex::Real denuc_weak_dT;
    amrex::Real denuc_weak_dYe;
    amrex::Real drate_Ni56_n_n_to_Ni58_approx_dYN{};
    amrex::Real drate_Ni58_to_Ni56_n_n_approx_dYN{};
    amrex::Real drate_Co53_n_n_to_Co55_approx_dYN{};
    amrex::Real drate_Co55_to_Co53_n_n_approx_dYN{};
    amrex::Real drate_Fe50_n_n_to_Fe52_approx_dYN{};
    amrex::Real drate_Fe52_to_Fe50_n_n_approx_dYN{};
    amrex::Real drate_Fe52_n_n_to_Fe54_approx_dYN{};
    amrex::Real drate_Fe54_to_Fe52_n_n_approx_dYN{};
    amrex::Real drate_Fe54_n_n_to_Fe56_approx_dYN{};
    amrex::Real drate_Fe56_to_Fe54_n_n_approx_dYN{};
    amrex::Real drate_Mn49_n_n_to_Mn51_approx_dYN{};
    amrex::Real drate_Mn51_to_Mn49_n_n_approx_dYN{};
    amrex::Real drate_Mn51_n_n_to_Mn53_approx_dYN{};
    amrex::Real drate_Mn53_to_Mn51_n_n_approx_dYN{};
    amrex::Real drate_Cr46_n_n_to_Cr48_approx_dYN{};
    amrex::Real drate_Cr48_to_Cr46_n_n_approx_dYN{};
    amrex::Real drate_Cr48_n_n_to_Cr50_approx_dYN{};
    amrex::Real drate_Cr50_to_Cr48_n_n_approx_dYN{};
    amrex::Real drate_Cr50_n_n_to_Cr52_approx_dYN{};
    amrex::Real drate_Cr52_to_Cr50_n_n_approx_dYN{};
    amrex::Real drate_V45_n_n_to_V47_approx_dYN{};
    amrex::Real drate_V47_to_V45_n_n_approx_dYN{};
    amrex::Real drate_V47_n_n_to_V49_approx_dYN{};
    amrex::Real drate_V49_to_V47_n_n_approx_dYN{};
    amrex::Real drate_Ti42_n_n_to_Ti44_approx_dYN{};
    amrex::Real drate_Ti44_to_Ti42_n_n_approx_dYN{};
    amrex::Real drate_Ti44_n_n_to_Ti46_approx_dYN{};
    amrex::Real drate_Ti46_to_Ti44_n_n_approx_dYN{};
    amrex::Real drate_Sc41_n_n_to_Sc43_approx_dYN{};
    amrex::Real drate_Sc43_to_Sc41_n_n_approx_dYN{};
    amrex::Real drate_Ca38_n_n_to_Ca40_approx_dYN{};
    amrex::Real drate_Ca40_to_Ca38_n_n_approx_dYN{};
    amrex::Real drate_Ca40_n_n_to_Ca42_approx_dYN{};
    amrex::Real drate_Ca42_to_Ca40_n_n_approx_dYN{};
    amrex::Real drate_Ar34_n_n_to_Ar36_approx_dYN{};
    amrex::Real drate_Ar36_to_Ar34_n_n_approx_dYN{};
    amrex::Real drate_Ar36_n_n_to_Ar38_approx_dYN{};
    amrex::Real drate_Ar38_to_Ar36_n_n_approx_dYN{};
    amrex::Real drate_Si28_n_n_to_Si30_approx_dYN{};
    amrex::Real drate_Si30_to_Si28_n_n_approx_dYN{};
    amrex::Real drate_S30_n_n_to_S32_approx_dYN{};
    amrex::Real drate_S32_to_S30_n_n_approx_dYN{};
    amrex::Real drate_S32_n_n_to_S34_approx_dYN{};
    amrex::Real drate_S34_to_S32_n_n_approx_dYN{};
};


#endif