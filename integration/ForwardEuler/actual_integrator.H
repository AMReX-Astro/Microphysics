#ifndef actual_integrator_H
#define actual_integrator_H

#include <network.H>
#include <actual_network.H>
#include <actual_rhs.H>
#include <burn_type.H>
#include <rate_type.H>
#include <temperature_integration.H>
#include <eos_type.H>
#include <eos.H>
#include <extern_parameters.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real calculate_dt (burn_t& state, Array1D<Real, 1, neqs>& ydot)
{
    // Our timestepping strategy is to prevent any quantity
    // from changing by more than a certain factor in any
    // timestep. We ignore this for species below atol_spec.

    Real dt = 1.0e200_rt;

    for (int n = 1; n <= NumSpec; ++n) {

        if (state.xn[n-1] >= atol_spec) {

            Real target_dX;
            if (ydot(n) > 0.0) {
                target_dX = (maximum_timestep_change_factor - 1.0_rt) * state.xn[n-1];
            } else {
                target_dX = (1.0_rt - 1.0_rt / maximum_timestep_change_factor) * state.xn[n-1];
            }

            Real dXdt = amrex::max(std::abs(ydot(n)), 1.0e-30_rt);

            dt = amrex::min(dt, target_dX / dXdt);

        }

    }

    if (integrate_temperature) {

        Real target_dT;
        if (ydot(net_itemp) > 0.0) {
            target_dT = (maximum_timestep_change_factor - 1.0_rt) * state.T;
        } else {
            target_dT = (1.0_rt - 1.0_rt / maximum_timestep_change_factor) * state.T;
        }

        Real dTdt = amrex::max(std::abs(ydot(net_itemp)), 1.0e-30_rt);

        dt = amrex::min(dt, target_dT / dTdt);

    }

    if (integrate_energy) {

        Real target_de;
        if (ydot(net_ienuc) > 0.0) {
            target_de = (maximum_timestep_change_factor - 1.0_rt) * state.e;
        } else {
            target_de = (1.0_rt - 1.0_rt / maximum_timestep_change_factor) * state.e;
        }

        Real dedt = amrex::max(std::abs(ydot(net_ienuc)), 1.0e-30_rt);

        dt = amrex::min(dt, target_de / dedt);

    }

    dt = amrex::min(dt, ode_max_dt);

    return dt;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void clean_state (burn_t& state)
{
    // Renormalize the abundances.

    normalize_abundances_burn(state);

    // Ensure that the temperature always stays within reasonable limits.

    state.T = amrex::min(MAX_TEMP, amrex::max(state.T, EOSData::mintemp));

}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void actual_integrator (burn_t& state, Real dt)
{
    // We assume that (rho, T) coming in are valid, do an EOS call
    // to fill the rest of the thermodynamic variables.

    eos_t eos_state;

    burn_to_eos(state, eos_state);

    eos(eos_input_rt, eos_state);

    eos_to_burn(eos_state, state);

    clean_state(state);

    state.self_heat = true;

    state.success = true;

    state.n_rhs = 0;

    Real t = 0.0;

    // When checking the integration time to see if we're done,
    // be careful with roundoff issues.

    const Real timestep_safety_factor = 1.0e-12_rt;

    int num_timesteps = 0;

    while (t < (1.0_rt - timestep_safety_factor) * dt && num_timesteps < ode_max_steps) {

        // Evaluate the RHS.

        Array1D<Real, 1, neqs> ydot;

        actual_rhs(state, ydot);

        state.n_rhs += 1;

        // Scale species terms by A (they come from the RHS in terms of Y, not X).

        for (int n = 1; n <= NumSpec; ++n) {
            ydot(n) *= aion[n-1];
        }

        // Calculate the timestep.

        Real dt_sub = calculate_dt(state, ydot);

        // Prevent the timestep from overshooting the final time.

        if (t + dt_sub > dt) {
            dt_sub = dt - t;
        }

        for (int n = 1; n <= NumSpec; ++n) {
            state.xn[n-1] += ydot(n) * dt_sub;
        }

        if (integrate_temperature) {
            state.T += ydot(net_itemp) * dt_sub;
        }
        if (integrate_energy) {
            state.e += ydot(net_ienuc) * dt_sub;
        }

        clean_state(state);

        t += dt_sub;
        ++num_timesteps;

    }

    if (num_timesteps >= ode_max_steps) {
        state.success = false;
    }

#ifndef AMREX_USE_CUDA
    if (burner_verbose) {
        // Print out some integration statistics, if desired.
        std::cout <<  "integration summary: " << std::endl;
        std::cout <<  "dens: " << state.rho << " temp: " << state.T << std::endl;
        std::cout <<  "energy released: " << state.e - eos_state.e << std::endl;
        std::cout <<  "number of steps taken: " << num_timesteps << std::endl;
        std::cout <<  "number of f evaluations: " << state.n_rhs << std::endl;
    }
#endif

    // If we failed, print out the current state of the integration.

    if (!state.success) {
#ifndef AMREX_USE_CUDA
        std::cout << "ERROR: integration failed in net" << std::endl;
        std::cout << "time = " << t << std::endl;
        std::cout << "dens = " << state.rho << std::endl;
        std::cout << "temp start = " << eos_state.T << std::endl;
        std::cout << "xn start = ";
        for (int n = 0; n < NumSpec; ++n) {
            std::cout << eos_state.xn[n] << " ";
        }
        std::cout << std::endl;
        std::cout << "temp current = " << state.T << std::endl;
        std::cout << "xn current = ";
        for (int n = 0; n < NumSpec; ++n) {
            std::cout << state.xn[n] << " ";
        }
        std::cout << std::endl;
        std::cout << "energy generated = " << state.e - eos_state.e << std::endl;
#endif
    }

}

#endif
