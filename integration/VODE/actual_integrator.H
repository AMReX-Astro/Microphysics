#ifndef actual_integrator_H
#define actual_integrator_H

// Common variables and routines for burners
// that use VODE for their integration.

#include <network.H>
#include <burn_type.H>
#include <eos_type.H>
#include <eos.H>
#include <extern_parameters.H>
#include <vode_type.H>
#include <vode_dvode.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void actual_integrator (burn_t& state, Real dt, bool is_retry=false)
{
    dvode_t vode_state;

    // Set the tolerances.

    if (!is_retry) {
        for (int n = 1; n <= NumSpec; ++n) {
            vode_state.atol(n) = atol_spec; // mass fractions
        }
        vode_state.atol(net_ienuc) = atol_enuc; // energy generated

        for (int n = 1; n <= NumSpec; ++n) {
            vode_state.rtol(n) = rtol_spec; // mass fractions
        }
        vode_state.rtol(net_ienuc) = rtol_enuc; // energy generated
    } else {
        for (int n = 1; n <= NumSpec; ++n) {
            vode_state.atol(n) = retry_atol_spec; // mass fractions
        }
        vode_state.atol(net_ienuc) = retry_atol_enuc; // energy generated

        for (int n = 1; n <= NumSpec; ++n) {
            vode_state.rtol(n) = retry_rtol_spec; // mass fractions
        }
        vode_state.rtol(net_ienuc) = retry_rtol_enuc; // energy generated
    }

    // set the Jacobian type
    if (!(is_retry && retry_swap_jacobian)) {
        vode_state.jacobian_type = jacobian;
    } else {
        vode_state.jacobian_type = (jacobian == 1) ? 2 : 1;
    }

    // Start off by assuming a successful burn.

    state.success = true;

    // Initialize the integration time.

    vode_state.t = 0.0_rt;
    vode_state.tout = dt;

    // Set the (inverse of the) timestep limiter.

    vode_state.HMXI = 1.0_rt / ode_max_dt;

    // We assume that (rho, T) coming in are valid, do an EOS call
    // to fill the rest of the thermodynamic variables.

    eos(eos_input_rt, state);

    // Fill in the initial integration state.

    burn_to_vode(state, vode_state);

    // Save the initial composition, temperature, and energy for our later diagnostics.

#ifndef AMREX_USE_GPU
    Real xn_in[NumSpec];
    for (int n = 0; n < NumSpec; ++n) {
        xn_in[n] = state.xn[n];
    }
    Real T_in = state.T;
#endif
    Real e_in = state.e;

    state.self_heat = true;

    // Call the integration routine.

    int istate = dvode(state, vode_state);

    // Copy the integration data back to the burn state.

    vode_to_burn(vode_state, state);

#ifdef NSE
    // compute the temperature based on the energy release -- we need 
    // this in case we failed in our burn here because we entered NSE

    eos(eos_input_re, state);

#endif

    // Subtract off the initial energy (the application codes expect
    // to get back only the generated energy during the burn).

    state.e -= e_in;

    // Normalize the final abundances.

    normalize_abundances_burn(state);

    // Get the number of RHS and Jacobian evaluations.

    state.n_rhs = vode_state.NFE;
    state.n_jac = vode_state.NJE;

    // VODE does not always fail even though it can lead to unphysical states.
    // Add some checks that indicate a burn fail even if VODE thinks the
    // integration was successful.

    if (istate < 0) {
        state.success = false;
    }

    for (int n = 1; n <= NumSpec; ++n) {
        if (vode_state.y(n) < -vode_failure_tolerance) {
            state.success = false;
        }

        if (vode_state.y(n) > 1.0_rt + vode_failure_tolerance) {
            state.success = false;
        }
    }

#ifndef AMREX_USE_GPU
    if (burner_verbose) {
        // Print out some integration statistics, if desired.
        std::cout <<  "integration summary: " << std::endl;
        std::cout <<  "dens: " << state.rho << " temp: " << state.T << std::endl;
        std::cout << " energy released: " << state.e << std::endl;
        std::cout <<  "number of steps taken: " << vode_state.NST << std::endl;
        std::cout <<  "number of f evaluations: " << vode_state.NFE << std::endl;
    }
#endif

    // If we failed, print out the current state of the integration.

    if (!state.success) {
        if (istate != -100) {
#ifndef AMREX_USE_GPU
            std::cout << Font::Bold << FGColor::Red << "[ERROR] integration failed in net" << ResetDisplay << std::endl;
            std::cout << "istate = " << istate << std::endl;
            if (istate == 2) {
                std::cout << "  VODE exited successfully, but a check on the data values failed" << std::endl;
            }
            std::cout << "time = " << vode_state.t << std::endl;
            std::cout << "dens = " << state.rho << std::endl;
            std::cout << "temp start = " << T_in << std::endl;
            std::cout << "xn start = ";
            for (int n = 0; n < NumSpec; ++n) {
                std::cout << xn_in[n] << " ";
            }
            std::cout << std::endl;
            std::cout << "temp current = " << state.T << std::endl;
            std::cout << "xn current = ";
            for (int n = 0; n < NumSpec; ++n) {
                std::cout << state.xn[n] << " ";
            }
            std::cout << std::endl;
            std::cout << "energy generated = " << state.e << std::endl;
#endif
        } else {
#ifndef AMREX_USE_GPU
            std::cout << "burn entered NSE during integration" << std::endl;
#endif
        }
    }

}

#endif
