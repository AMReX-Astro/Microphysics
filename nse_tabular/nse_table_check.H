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
#include <eos.H>


AMREX_GPU_HOST_DEVICE AMREX_INLINE
bool in_nse(burn_t& state, const bool relax = false) {

    using namespace Species;

    Real f = relax ? nse_relax_factor : 1.0_rt;

    state.nse = false;

    if (state.rho > f * rho_nse && state.T > f * T_nse) {

        // Ma checks on Fe-group (for our composition, this means Cr48, Fe52, Fe54, Ni56)
        // and C-group (for us, that is C12, N14)
        // and He-group (for us, that is H1, He3, He4)
        // also make sure there is not a lot of O16 (O16 burning will dominate then) or Si28

#ifdef STRANG
        Real Fe_group = state.xn[Cr48-1] + state.xn[Fe52-1] +
                        state.xn[Fe54-1] + state.xn[Ni56-1];
        Real C_group = state.xn[C12-1] + state.xn[N14-1];
        Real He_group = state.xn[H1-1] + state.xn[He3-1] + state.xn[He4-1];
        Real O_group = state.xn[O16-1];
        Real Si_group = state.xn[Si28-1];
#else
        // we need to get the mass fractions from the conserved state passed in

        Real Fe_group = (state.y[SFS+Cr48-1] + state.y[SFS+Fe52-1] +
                         state.y[SFS+Fe54-1] + state.y[SFS+Ni56-1]) / state.rho;
        Real C_group = (state.y[SFS+C12-1] + state.y[SFS+N14-1]) / state.rho;
        Real He_group = (state.y[SFS+H1-1] + state.y[SFS+He3-1] + state.y[SFS+He4-1]) / state.rho;
        Real O_group = state.y[SFS+O16-1] / state.rho;
        Real Si_group = state.y[SFS+Si28-1] / state.rho;
#endif

        if (Fe_group + He_group > f * He_Fe_nse &&
            C_group < C_nse / f &&
            O_group < O_nse / f &&
            Si_group < Si_nse / f) {
            state.nse = true;
        }
    }

    return state.nse;

}

#endif
