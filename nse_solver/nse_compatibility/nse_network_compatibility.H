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
    burn_state.e = eos_state.e;

    burn_state.y[SRHO] = eos_state.rho;
    burn_state.y[SMX] = 0.0;
    burn_state.y[SMY] = 0.0;
    burn_state.y[SMZ] = 0.0;

    burn_state.y[SEINT] = eos_state.rho * eos_state.e;
    burn_state.y[SEDEN] = burn_state.y[SEINT];

    burn_state.ydot_a[SRHO] = 0.0;
    burn_state.ydot_a[SEINT] = 0.0;

    for (int n = 0; n < NumSpec; ++n) {
        burn_state.xn[n] = eos_state.xn[n];
        burn_state.y[SFS+n] = eos_state.xn[n] * burn_state.rho;
    }

    burn_state.y_e = eos_state.y_e;

    burn_state.sdc_iter = 1;
    burn_state.num_sdc_iters = 1;

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
    amrex::Real dlogt = 0.0_rt;
    if (nsteps == 1) {
        dlogt = (std::log10(tmax) - std::log10(tfirst));
    } else {
        dlogt = (std::log10(tmax) - std::log10(tfirst)) / (nsteps - 1);
    }

    // save the initial state -- we'll use this to determine
    // how much things changed over the entire burn

    burn_t burn_state_in = burn_state;

    amrex::Real t = 0.0;

    // store the initial internal energy -- we'll update this after
    // each substep

    amrex::Real energy_initial = burn_state.e;

    // Now since we compiled with USE_NSE_NET=TRUE,
    // it will try to burn using NSE, to prevent that, set extremely small tolerance

    ase_tol = 1.e-100;
    burn_state.mu_p = mu_p;
    burn_state.mu_n = mu_n;

    // loop over steps, burn, and output the current state

    for (int n = 0; n < nsteps; n++){

        // compute the time we wish to integrate to

        amrex::Real tend = std::pow(10.0_rt, std::log10(tfirst) + dlogt * n);
        amrex::Real dt = tend - t;

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
    auto nse_state = get_actual_nse_state(nse_input_rty, burn_state, 1.e-12);

    // Now print out the result.
    std::cout << std::scientific << std::setprecision(8);
    std::cout << "------------------------------------" << std::endl;
    std::cout << " - added e = " << burn_state.e - burn_state_in.e << std::endl;
    std::cout << " Initial vs. Final T = " << burn_state_in.T << "  "
              << burn_state.T << std::endl;
    std::cout << " Initial vs. Final density = " << burn_state_in.rho << "  "
              << burn_state.rho << std::endl;
    std::cout << " Initial vs. Final vs. NSE y_e = " << burn_state_in.y_e << "  "
              << burn_state.y_e << "  " <<  nse_state.y_e << std::endl;
    std::cout << " Solved Chemical Potentials are: " << burn_state.mu_p << " and "
              << burn_state.mu_n << std::endl;
    std::cout << "------------------------------------" << std::endl;
    std::cout << "Element" << std::setw(14) << "Burn Xs" << std::setw(20) << "NSE Xs "
              << std::setw(20) << "abs err" << std::setw(20) << "rel err" << std::endl;
    for (int n = 0; n < NumSpec; ++n) {
        const std::string& element = short_spec_names_cxx[n];
        amrex::Real abs_err = std::abs(burn_state.xn[n] - nse_state.xn[n]);
        amrex::Real rel_err = abs_err / burn_state.xn[n];
        std::cout << "  " << element << std::setw(20-element.size())
                  << burn_state.xn[n] << std::setw(20)
                  << nse_state.xn[n] << std::setw(20) << abs_err << std::setw(20)
                  << rel_err << std::endl;
    }
    std::cout << "------------------------------------" << std::endl;
}

void nse_network_compatibility()
{
    // This function compares the integrated equilibrium mass fractions
    // with the results from NSE calculations given any network.
    // Note that this network should not have any weak rates.

    eos_t eos_state;
    eos_state.T = T;
    eos_state.rho = rho;
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = 1.0_rt;
    }
    normalize_abundances_burn(eos_state);
    eos(eos_input_rt, eos_state);
    burn_cell_c(eos_state);
}
#endif
