#ifndef actual_integrator_H
#define actual_integrator_H

#include <iomanip>

#include <network.H>
#include <burn_type.H>
#include <eos_type.H>
#include <eos.H>
#include <extern_parameters.H>
#include <rock4_type.H>
#include <rock4.H>

template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void actual_integrator (BurnT& state, Real dt)
{
    rock4_t rock_state{};

    // Set the tolerances.

    rock_state.atol_spec = atol_spec; // mass fractions
    rock_state.atol_enuc = atol_enuc; // energy generated

    rock_state.rtol_spec = rtol_spec; // mass fractions
    rock_state.rtol_enuc = rtol_enuc; // energy generated

    // Start off by assuming a successful burn.

    state.success = true;

    // Initialize the integration time.

    rock_state.t = 0.0_rt;
    rock_state.tend = dt;

    // We assume that (rho, T) coming in are valid, do an EOS call
    // to fill the rest of the thermodynamic variables.

    eos(eos_input_rt, state);

    // Fill in the initial integration state.

    burn_to_integrator(state, rock_state);

    // Save the initial composition, temperature, and energy for our later diagnostics.

#ifndef AMREX_USE_GPU
    Real xn_in[NumSpec];
    for (int n = 0; n < NumSpec; ++n) {
        xn_in[n] = state.xn[n];
    }
    Real T_in = state.T;
#endif
    Real e_in = state.e;

    // Call the integration routine.

    int ierr = rock4(state, rock_state);

    // Copy the integration data back to the burn state.

    integrator_to_burn(rock_state, state);

    // Subtract off the initial energy (the application codes expect
    // to get back only the generated energy during the burn).
    // Don't subtract it for primordial chem
    if (subtract_internal_energy) {
        state.e -= e_in;
    }

    // Normalize the final abundances.
    // Don't normalize for primordial chem
    if (!use_number_densities) {
        normalize_abundances_burn(state);
    }

    // Get the number of RHS and Jacobian evaluations.

    state.n_rhs = rock_state.nfe;
    state.n_jac = 0;
    state.n_step = rock_state.nsteps;

    if (ierr > 2) {
        state.success = false;
    }

    // ROCK4 does not always fail even though it can lead to unphysical states.
    // Add some checks that indicate a burn fail even if ROCK4 thinks the
    // integration was successful.

    for (int n = 1; n <= NumSpec; ++n) {
        if (rock_state.y(n) < -rock_failure_tolerance) {
            state.success = false;
        }

        // Don't enforce the condition below
        // for primordial chem
        if (!use_number_densities) {
            if (rock_state.y(n) > 1.0_rt + rock_failure_tolerance) {
                state.success = false;
            }
        }
    }

#ifndef AMREX_USE_GPU
    if (burner_verbose) {
        // Print out some integration statistics, if desired.
        std::cout <<  "integration summary: " << std::endl;
        std::cout <<  "dens: " << state.rho << " temp: " << state.T << std::endl;
        std::cout << " energy released: " << state.e << std::endl;
        std::cout <<  "number of steps taken: " << rock_state.nsteps << std::endl;
        std::cout <<  "number of f evaluations: " << rock_state.nfe << std::endl;
    }
#endif

    // If we failed, print out the current state of the integration.

    if (!state.success) {
#ifndef AMREX_USE_GPU
        std::cout << Font::Bold << FGColor::Red << "[ERROR] integration failed in net" << ResetDisplay << std::endl;
        std::cout << "ierr = " << ierr << std::endl;
        std::cout << "zone = (" << state.i << ", " << state.j << ", " << state.k << ")" << std::endl;
        std::cout << "time = " << rock_state.t << std::endl;
        std::cout << "dt = " << std::setprecision(16) << dt << std::endl;
        std::cout << "temp start = " << std::setprecision(16) << T_in << std::endl;
        std::cout << "xn start = ";
        for (int n = 0; n < NumSpec; ++n) {
            std::cout << std::setprecision(16) << xn_in[n] << " ";
        }
        std::cout << std::endl;
        std::cout << "dens current = " << std::setprecision(16) << state.rho << std::endl;
        std::cout << "temp current = " << std::setprecision(16) << state.T << std::endl;
        std::cout << "xn current = ";
        for (int n = 0; n < NumSpec; ++n) {
            std::cout << std::setprecision(16) << state.xn[n] << " ";
        }
        std::cout << std::endl;
        std::cout << "energy generated = " << state.e << std::endl;
#endif
    }
}

#endif
