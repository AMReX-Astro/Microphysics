#ifndef _vode_integrator_H_
#define _vode_integrator_H_

// Common variables and routines for burners
// that use VODE for their integration.

#include <network.H>
#include <burn_type.H>
#include <temperature_integration.H>
#include <eos_type.H>
#include <eos.H>
#include <extern_parameters.H>
#include <vode_type.H>
#include <vode_dvode.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void vode_integrator (burn_t& state, Real dt)
{

    dvode_t vode_state;

    const Real failure_tolerance = 1.e-2_rt;

    // Set the tolerances.  We will be more relaxed on the temperature
    // since it is only used in evaluating the rates.
    //
    // **NOTE** if you reduce these tolerances, you probably will need
    // to (a) decrease dT_crit, (b) increase the maximum number of
    // steps allowed.

    for (int n = 1; n <= NumSpec; ++n) {
        vode_state.atol(n) = atol_spec; // mass fractions
    }
    vode_state.atol(net_itemp) = atol_temp; // temperature
    vode_state.atol(net_ienuc) = atol_enuc; // energy generated

    for (int n = 1; n <= NumSpec; ++n) {
        vode_state.rtol(n) = rtol_spec; // mass fractions
    }
    vode_state.rtol(net_itemp) = rtol_temp; // temperature
    vode_state.rtol(net_ienuc) = rtol_enuc; // energy generated

    // Start off by assuming a successful burn.

    state.success = true;

    // Initialize the integration time.

    vode_state.t = 0.0_rt;
    vode_state.tout = dt;

    // We assume that (rho, T) coming in are valid, do an EOS call
    // to fill the rest of the thermodynamic variables.

    eos_t eos_state;

    burn_to_eos(state, eos_state);

    eos(eos_input_rt, eos_state);

    eos_to_burn(eos_state, state);

    // Fill in the initial integration state.

    burn_to_vode(state, vode_state);

    // Save the initial energy for our later diagnostics.

    Real e_in = state.e;

    // If we are using the dT_crit functionality and therefore doing a linear
    // interpolation of the specific heat in between EOS calls, do a second
    // EOS call here to establish an initial slope.

    state.T_old = state.T;
    state.cv_old = state.cv;
    state.cp_old = state.cp;

    if (dT_crit < 1.0e19_rt) {

        eos_state.T *= (1.0_rt + std::sqrt(std::numeric_limits<Real>::epsilon()));

        eos(eos_input_rt, eos_state);

        state.dcvdT = (eos_state.cv - state.cv_old) / (eos_state.T - state.T_old);
        state.dcpdT = (eos_state.cp - state.cp_old) / (eos_state.T - state.T_old);

    }

    if (burning_mode == 0) {
        state.self_heat = false;
    }
    else if (burning_mode == 1) {
        state.self_heat = true;
    }
#ifndef AMREX_USE_GPU
    else {
        amrex::Error("Unsupported burning_mode");
    }
#endif

    // Call the integration routine.

    int istate = dvode(state, vode_state);

    // VODE does not always fail even though it can lead to unphysical states,
    // so add some sanity checks that trigger a retry even if VODE thinks
    // the integration was successful.

    if (istate < 0) {
        state.success = false;
    }

    if (vode_state.y(net_itemp) < 0.0_rt) {
        state.success = false;
    }

    for (int n = 1; n <= NumSpec; ++n) {
        if (vode_state.y(n) < -failure_tolerance) {
            state.success = false;
        }

        if (vode_state.y(n) > 1.e0_rt + failure_tolerance) {
            state.success = false;
        }
    }


    // If we failed, print out the current state of the integration.

    if (!state.success) {
#ifndef AMREX_USE_CUDA
        std::cout << "ERROR: integration failed in net" << std::endl;
        std::cout << "istate = " << istate << std::endl;
        std::cout << "time = " << vode_state.t << std::endl;
#endif

        state.success = false;
        return;
    }

    // Subtract the energy offset

    vode_state.y(net_ienuc) -= e_in;

    // Copy the integration data back to the burn state.

    vode_to_burn(vode_state, state);

    // Get the number of RHS and Jacobian evaluations

    state.n_rhs = vode_state.NFE;
    state.n_jac = vode_state.NJE;

    normalize_abundances_burn(state);

#ifndef AMREX_USE_CUDA
    if (burner_verbose) {
        // Print out some integration statistics, if desired.
        std::cout <<  "integration summary: " << std::endl;
        std::cout <<  "dens: " << state.rho << " temp: " << state.T << std::endl;
        std::cout << " energy released: " << state.e << std::endl;
        std::cout <<  "number of steps taken: " << vode_state.NST << std::endl;
        std::cout <<  "number of f evaluations: " << vode_state.NFE << std::endl;
    }
#endif

}

#endif
