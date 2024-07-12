#ifndef BURN_CELL_H
#define BURN_CELL_H

#include <extern_parameters.H>
#include <eos.H>
#include <network.H>
#include <burner.H>
#include <fstream>
#include <iostream>
#include <nse_solver.H>
#include <cmath>

using namespace unit_test_rp;

AMREX_INLINE
void burn_cell_c()
{

    burn_t state;

    // set a reference cell size

    state.dx = 1.0e6_rt;

    amrex::Real dlogrho = (std::log10(rho_max) - std::log10(rho_min))/static_cast<amrex::Real>(nrho-1);
    amrex::Real dlogT = (std::log10(T_max) - std::log10(T_min))/static_cast<amrex::Real>(nT-1);
    amrex::Real dYe = (Ye_max - Ye_min)/(nye-1);

    for (int iye = 0; iye < nye; ++iye) {
        for (int irho = 0; irho < nrho; ++irho) {
            for (int itemp = 0; itemp < nT; ++itemp) {

                amrex::Real T = std::pow(10.0, std::log10(T_min) + itemp * dlogT);
                amrex::Real rho = std::pow(10.0, std::log10(rho_min) + irho * dlogrho);
                amrex::Real Ye = Ye_min + iye * dYe;
                state.T = T;
                state.rho = rho;
                state.y_e = Ye;

                if (state.y_e > 0.52_rt){
                    state.mu_p = -4.0_rt;
                    state.mu_n = -14.0_rt;
                }
                else if (state.y_e > 0.48_rt){
                    state.mu_p = -8.0_rt;
                    state.mu_n = -11.0_rt;
                }
                else if (state.y_e > 0.4_rt){
                    state.mu_p = -10.0_rt;
                    state.mu_n = -7.0_rt;
                }
                else{
                    state.mu_p = -15.0_rt;
                    state.mu_n = -5.0_rt;
                }

                // find the  nse state

                const bool assume_ye_is_valid = true;
                amrex::Real eps = 1.e-10_rt;

                use_hybrid_solver = 1;

                auto nse_state = get_actual_nse_state(state, eps, assume_ye_is_valid);

                for (int i = 0; i < NumSpec; ++i){
                  state.y[SFS+i] = nse_state.y[SFS+i];
                }

                bool in_nse_state = in_nse(state);

                std::cout << std::scientific;
                std::cout << std::setw(20) << state.rho << " "
                          << std::setw(20) << state.T << " " << std::fixed
                          << std::setw(20) << state.y_e << " "
                          << std::setw(20) << state.mu_p << " "
                          << std::setw(20) << state.mu_n << " "
                          << std::setw(20) << in_nse_state << std::endl;

            }
        }
    }
}
#endif
