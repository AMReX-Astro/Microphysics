
#include <extern_parameters.H>
#include <eos.H>
#include <network.H>
#include <burner.H>
#include <fstream>
#include <iostream>
#include <actual_nse.H>
#include <cmath>
#include <limits>

void burn_cell_c()
{

    Real min_Ye = std::numeric_limits<Real>::max();
    for (int n = 0; n < NumSpec; ++n) {
        min_Ye = amrex::min(min_Ye, zion[n]/aion[n]);
    }
    std::cout << "minimum Ye for this network = " << min_Ye << std::endl;

    burn_t state;

    Real dlogrho = (std::log10(rho_max) - std::log10(rho_min))/static_cast<Real>(nrho-1);
    Real dlogT = (std::log10(T_max) - std::log10(T_min))/static_cast<Real>(nT-1);
    Real dYe = (Ye_max - Ye_min)/(nye-1);

    for (int iye = 0; iye < nye; ++iye) {
        for (int irho = 0; irho < nrho; ++irho) {
            for (int itemp = 0; itemp < nT; ++itemp) {

                Real T = std::pow(10.0, std::log10(T_min) + itemp * dlogT);
                Real rho = std::pow(10.0, std::log10(rho_min) + irho * dlogrho);
                Real Ye = Ye_min + iye * dYe;

                state.T = T;
                state.rho = rho;
                state.y_e = Ye;

                // initial guess

                if (state.T >= 4.5e9_rt) {

                    if (state.y_e < 0.45_rt) {
                        state.mu_p = -15.0;
                        state.mu_n = -4.0;
                    } else if (state.y_e < 0.5_rt) {
                        state.mu_p = -10.0;
                        state.mu_n = -6.0;
                    } else {
                        state.mu_p = -3.0;
                        state.mu_n = -14.0;
                    }

                } else {
                    // lower temperature

                    if (state.y_e < 0.45_rt) {
                        if (state.rho < 1.e8_rt) {
                            state.mu_p = -14.0;
                            state.mu_n = -5.0;
                        } else {
                            state.mu_p = -16.3;
                            state.mu_n = -3.2;
                        }
                    } else if (state.y_e < 0.5_rt) {
                        if (state.rho < 1.e8_rt) {
                            state.mu_p = -6.0;
                            state.mu_n = -11.0;
                        } else {
                            state.mu_p = -9.9;
                            state.mu_n = -8.0;
                        }
                    } else {
                        state.mu_p = -4.0;
                        state.mu_n = -13.0;
                    }

                }

                // find the  nse state

                const bool assume_ye_is_valid = true;
                Real eps = 1.e-4;

                auto nse_state = get_actual_nse_state(state, eps, assume_ye_is_valid);
                std::cout << std::scientific;
                std::cout << std::setw(20) << state.rho << " "
                          << std::setw(20) << state.T << " " << std::fixed
                          << std::setw(20) << state.y_e << " "
                          << std::setw(20) << state.mu_p << " "
                          << std::setw(20) << state.mu_n << std::endl;

            }
        }
    }
}
