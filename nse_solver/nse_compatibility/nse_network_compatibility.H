#ifndef NSE_NETWORK_COMPATIBILITY_H
#define NSE_NETWORK_COMPATIBILITY_H

#include <extern_parameters.H>
#include <eos.H>
#include <network.H>
#include <burner.H>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <react_util.H>

using namespace unit_test_rp;

void burn_cell_c(const eos_t& eos_state)
{
    // Given eos_state with density, temperature and mass fractions

    burn_t burn_state;

    burn_state.rho = eos_state.rho;
    burn_state.T = eos_state.T;
    for (int n = 0; n < NumSpec; ++n) {
        burn_state.xn[n] = eos_state.xn[n];
    }
    burn_state.y_e = eos_state.y_e;

    burn_state.i = 0;
    burn_state.j = 0;
    burn_state.k = 0;
    burn_state.T_fixed = -1.0_rt;

    // normalize -- just in case
    normalize_abundances_burn(burn_state);

    // call the EOS to set initial e -- it actually doesn't matter to
    // the burn but we need to keep track of e to get a valid
    // temperature for the burn if we substep

    eos(eos_input_rt, burn_state);

    // we will divide the total integration time into nsteps that are
    // logarithmically spaced

    if (tfirst == 0.0_rt) {
        if (nsteps == 1) {
            tfirst = tmax;
        } else {
            tfirst = tmax / nsteps;
        }
    }
    Real dlogt = 0.0_rt;
    if (nsteps == 1) {
        dlogt = (std::log10(tmax) - std::log10(tfirst));
    } else {
        dlogt = (std::log10(tmax) - std::log10(tfirst)) / (nsteps - 1);
    }

    // save the initial state -- we'll use this to determine
    // how much things changed over the entire burn

    burn_t burn_state_in = burn_state;

    Real t = 0.0;

    // store the initial internal energy -- we'll update this after
    // each substep

    Real energy_initial = burn_state.e;

    // loop over steps, burn, and output the current state

    for (int n = 0; n < nsteps; n++){

        // compute the time we wish to integrate to

        Real tend = std::pow(10.0_rt, std::log10(tfirst) + dlogt * n);
        Real dt = tend - t;

        burner(burn_state, dt);

        if (! burn_state.success) {
            amrex::Error("integration failed");
        }

        // state.e represents the change in energy over the burn (for
        // just this sybcycle), so turn it back into a physical energy

        burn_state.e += energy_initial;

        // reset the initial energy for the next subcycle

        energy_initial = burn_state.e;

        t += dt;
    }

    // Now find the mass fraction via NSE calculation
    burn_state.mu_p = -3.0_rt;
    burn_state.mu_n = -12.0_rt;
    auto nse_state = get_actual_nse_state(burn_state);

    if (nse_dx_independent != 1 || (nse_dx_independent == 1 && !burn_state.nse)) {
        std::cout << std::scientific << std::setprecision(3);
        std::cout << "------------------------------------" << std::endl;
        std::cout << " - added e = " << burn_state.e - burn_state_in.e << std::endl;
        std::cout << " Initial vs. Final T = " << burn_state_in.T << "  "
                  << burn_state.T << std::endl;
        std::cout << " Initial vs. Final density = " << burn_state_in.rho << "  "
                  << burn_state.rho << std::endl;
        std::cout << " Initial vs. Final vs. NSE y_e = " << burn_state_in.y_e << "  "
                  << burn_state.y_e << "  " <<  nse_state.y_e << std::endl;
        std::cout << "------------------------------------" << std::endl;
        std::cout << "Element" << std::setw(14) << "Burn Xs" << std::setw(20) << "NSE Xs "
                  << std::setw(20) << "abs err" << std::setw(20) << "rel err" << std::endl;
        for (int n = 0; n < NumSpec; ++n) {
            const std::string& element = short_spec_names_cxx[n];
            Real abs_err = std::abs(burn_state.xn[n] - nse_state.xn[n]);
            Real rel_err = abs_err / burn_state.xn[n];
            std::cout << "  " << element << std::setw(20-element.size())
                      << burn_state.xn[n] << std::setw(20)
                      << nse_state.xn[n] << std::setw(20) << abs_err << std::setw(20)
                      << rel_err << std::endl;
        }
        std::cout << "------------------------------------" << std::endl;
    }
}

void nse_network_compatibility()
{
    // This function compares the integrated equilibrium mass fractions
    // with the results from NSE calculations given any network.
    // Note that this network should not have any weak rates.

    Real dlogrho = (std::log10(rho_max) - std::log10(rho_min))/static_cast<Real>(nrho-1);
    Real dlogT = (std::log10(T_max) - std::log10(T_min))/static_cast<Real>(nT-1);

    // Create massfraction 2d array to hold 3 differents sets of initial mass fractions
    // Corresponding to ye < 0.5, ye = 0.5, and ye > 0.5

    amrex::Array2D<Real, 1, 3, 1, NumSpec+1, Order::C> massfractions;

    bool condition1_met = false;
    bool condition2_met = false;
    bool condition3_met = false;

    for (int n = 0; n < NumSpec; ++n) {
        if (aion[n] == 2 * zion[n] && !condition1_met) {

            // Set the first encountered nuclei with Z/A = 0.5
            massfractions(1, n+1) = 0.7_rt;
            massfractions(2, n+1) = 1.0_rt;
            massfractions(3, n+1) = 0.7_rt;
            condition1_met = true;
            continue;
        }
        if (aion[n] > 2 * zion[n] && !condition2_met) {

            // This is for achieving ye < 0.5
            massfractions(1, n+1) = 0.3_rt;
            condition2_met = true;
            continue;
        }
        if (aion[n] < 2 * zion[n] && !condition3_met) {

            // This is for achieving ye > 0.5
            massfractions(3, n+1) = 0.3_rt;
            condition3_met = true;
            continue;
        }
        massfractions(1, n+1) = 0.0_rt;
        massfractions(2, n+1) = 0.0_rt;
        massfractions(3, n+1) = 0.0_rt;
    }

    eos_t eos_state;

    for (int iye = 0; iye < 3; ++iye) {
        for (int irho = 0; irho < nrho; ++irho) {
            for (int itemp = 0; itemp < nT; ++itemp) {

                Real T = std::pow(10.0, std::log10(T_min) + itemp * dlogT);
                Real rho = std::pow(10.0, std::log10(rho_min) + irho * dlogrho);

                eos_state.T = T;
                eos_state.rho = rho;
                for (int n = 0; n < NumSpec; ++n) {
                    eos_state.xn[n] = massfractions(iye+1, n+1);
                }
                eos(eos_input_rt, eos_state);
                burn_cell_c(eos_state);
            }
        }
    }
}
#endif
