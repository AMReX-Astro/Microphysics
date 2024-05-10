#ifndef NSE_TABLE_CHECK_H
#define NSE_TABLE_CHECK_H

#include <iostream>
#include <fstream>
#include <actual_network.H>

#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Array.H>
#include <AMReX_REAL.H>

#include <extern_parameters.H>

#include <burn_type.H>
#include <eos_type.H>


template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
bool in_nse(T& state, const bool relax = false) {

    using namespace Species;

    Real f = relax ? network_rp::nse_relax_factor : 1.0_rt;

    bool nse_check = false;

    if (state.rho > f * network_rp::rho_nse && state.T > f * network_rp::T_nse) {

        if (state.T > network_rp::T_always_in_nse) {
           nse_check = true;

        } else {

            // consider composition

            // Ma et al. 2013 checks on Fe-group (for our composition,
            // this means Cr48, Fe52, Fe54, Ni56)
            // and C-group (for us, that is C12, N14)
            // and He-group (for us, that is H1, He3, He4)

            // also make sure there is not a lot of O16 (O16 burning
            // will dominate then) or Si28

            Real Fe_group{};
            Real C_group{};
            Real He_group{};
            Real O_group{};
            Real Si_group{};

            if constexpr (std::is_same_v<T, burn_t>) {
                // for a burn_t, we need to use a different field
                // depending on whether we are strang of simplified-SDC

#ifdef STRANG
                Fe_group = state.xn[Cr48-1] + state.xn[Fe52-1] +
                           state.xn[Fe54-1] + state.xn[Ni56-1];
                C_group = state.xn[C12-1] + state.xn[N14-1];
                He_group = state.xn[H1-1] + state.xn[He3-1] + state.xn[He4-1];
                O_group = state.xn[O16-1];
                Si_group = state.xn[Si28-1];
#else
                // we need to get the mass fractions from the conserved state passed in

                Fe_group = (state.y[SFS+Cr48-1] + state.y[SFS+Fe52-1] +
                            state.y[SFS+Fe54-1] + state.y[SFS+Ni56-1]) / state.rho;
                C_group = (state.y[SFS+C12-1] + state.y[SFS+N14-1]) / state.rho;
                He_group = (state.y[SFS+H1-1] + state.y[SFS+He3-1] + state.y[SFS+He4-1]) / state.rho;
                O_group = state.y[SFS+O16-1] / state.rho;
                Si_group = state.y[SFS+Si28-1] / state.rho;
#endif

            } else {

                // this covers the various eos_t's
                Fe_group = state.xn[Cr48-1] + state.xn[Fe52-1] +
                           state.xn[Fe54-1] + state.xn[Ni56-1];
                C_group = state.xn[C12-1] + state.xn[N14-1];
                He_group = state.xn[H1-1] + state.xn[He3-1] + state.xn[He4-1];
                O_group = state.xn[O16-1];
                Si_group = state.xn[Si28-1];
            }

            if (Fe_group + He_group > f * network_rp::He_Fe_nse &&
                C_group < network_rp::C_nse / f &&
                O_group < network_rp::O_nse / f &&
                Si_group < network_rp::Si_nse / f) {

                nse_check = true;
            }

        }

    }

    if constexpr (std::is_same<T, burn_t>::value) {
        state.nse = nse_check;
    }

    return nse_check;

}

#endif
