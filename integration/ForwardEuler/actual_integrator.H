#ifndef actual_integrator_H
#define actual_integrator_H

#include <network.H>
#include <actual_network.H>
#include <rhs.H>
#include <actual_rhs.H>
#include <burn_type.H>
#include <rate_type.H>
#include <eos_type.H>
#include <eos.H>
#include <extern_parameters.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real calculate_dt (burn_t& state, Array1D<Real, 1, NumSpec>& spec_rhs, Real& ener_rhs)
{
    // Our timestepping strategy is to prevent any quantity
    // from changing by more than a certain factor in any
    // timestep. We ignore this for species below atol_spec.

    Real dt = 1.0e200_rt;

    for (int n = 1; n <= NumSpec; ++n) {

        if (state.xn[n-1] >= atol_spec) {

            Real target_dX;
            if (spec_rhs(n) > 0.0) {
                target_dX = (maximum_timestep_change_factor - 1.0_rt) * state.xn[n-1];
            } else {
                target_dX = (1.0_rt - 1.0_rt / maximum_timestep_change_factor) * state.xn[n-1];
            }

            Real dXdt = amrex::max(std::abs(spec_rhs(n)), 1.0e-30_rt);

            dt = amrex::min(dt, target_dX / dXdt);

        }

    }

    if (integrate_energy) {

        Real target_de;
        if (ener_rhs > 0.0) {
            target_de = (maximum_timestep_change_factor - 1.0_rt) * state.e;
        } else {
            target_de = (1.0_rt - 1.0_rt / maximum_timestep_change_factor) * state.e;
        }

        Real dedt = amrex::max(std::abs(ener_rhs), 1.0e-30_rt);

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

    // Evaluate the EOS to get T from e.

    if (call_eos_in_rhs) {
        eos(eos_input_re, state);
    }

    // Ensure that the temperature always stays within reasonable limits.

    state.T = amrex::min(MAX_TEMP, amrex::max(state.T, EOSData::mintemp));

}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void evaluate_rhs (burn_t& state, Array1D<Real, 1, NumSpec>& spec_rhs, Real& ener_rhs)
{
    // Evaluate the RHS.

    Array1D<Real, 1, neqs> ydot;

    actual_rhs(state, ydot);

    state.n_rhs += 1;

    // Scale species terms by A (they come from the RHS in terms of Y, not X).

    for (int n = 1; n <= NumSpec; ++n) {
        spec_rhs(n) = aion[n-1] * ydot(n);
    }

    ener_rhs = ydot(net_ienuc);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void initialize_state (burn_t& state)
{
    // We assume that (rho, T) coming in are valid, do an EOS call
    // to fill the rest of the thermodynamic variables.

    eos(eos_input_rt, state);

    state.success = true;

    state.n_rhs = 0;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void actual_integrator (burn_t& state, Real dt)
{
    initialize_state(state);

    Real T_in = state.T;
    Real e_in = state.e;
    Real xn_in[NumSpec];
    for (int n = 0; n < NumSpec; ++n) {
        xn_in[n] = state.xn[n];
    }

    Real t = 0.0;

    // When checking the integration time to see if we're done,
    // be careful with roundoff issues.

    const Real timestep_safety_factor = 1.0e-12_rt;

    int num_timesteps = 0;

    while (t < (1.0_rt - timestep_safety_factor) * dt && num_timesteps < ode_max_steps) {

        // Evaluate the RHS.

        Array1D<Real, 1, NumSpec> spec_rhs;
        Real ener_rhs;

        evaluate_rhs(state, spec_rhs, ener_rhs);

       // Calculate the timestep.

        Real dt_sub = calculate_dt(state, spec_rhs, ener_rhs);

        // Prevent the timestep from overshooting the final time.

        if (t + dt_sub > dt) {
            dt_sub = dt - t;
        }

        // Apply the update.

        for (int n = 1; n <= NumSpec; ++n) {
            state.xn[n-1] += spec_rhs(n) * dt_sub;
        }

        if (integrate_energy) {
            state.e += ener_rhs * dt_sub;
        }

        clean_state(state);

        t += dt_sub;
        ++num_timesteps;

    }

    if (num_timesteps >= ode_max_steps) {
        state.success = false;
    }

    state.n_step = num_timesteps;
    state.time = t;

#ifndef AMREX_USE_GPU
    if (burner_verbose) {
        // Print out some integration statistics, if desired.
        std::cout <<  "integration summary: " << std::endl;
        std::cout <<  "dens: " << state.rho << " temp: " << state.T << std::endl;
        std::cout <<  "energy released: " << state.e - e_in << std::endl;
        std::cout <<  "number of steps taken: " << num_timesteps << std::endl;
        std::cout <<  "number of f evaluations: " << state.n_rhs << std::endl;
    }
#endif

    // If we failed, print out the current state of the integration.

    if (!state.success) {
#ifndef AMREX_USE_GPU
        std::cout << "ERROR: integration failed in net" << std::endl;
        std::cout << "time = " << t << std::endl;
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
        std::cout << "energy generated = " << state.e - e_in << std::endl;
#endif
    }

}

#endif
