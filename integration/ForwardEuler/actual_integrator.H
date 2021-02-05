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
Real calculate_dt (burn_t& state, Array1D<Real, 1, NumSpec>& dXdt, Real& dedt)
{
    // Our timestepping strategy is to prevent any quantity
    // from changing by more than a certain factor in any
    // timestep. We ignore this for species below atol_spec.

    Real dt = 1.0e200_rt;

    for (int n = 1; n <= NumSpec; ++n) {

        if (state.xn[n-1] >= atol_spec) {

            Real target_dX;
            if (dXdt(n) > 0.0) {
                target_dX = (maximum_timestep_change_factor - 1.0_rt) * state.xn[n-1];
            } else {
                target_dX = (1.0_rt - 1.0_rt / maximum_timestep_change_factor) * state.xn[n-1];
            }

            Real dXdt_abs = amrex::max(std::abs(dXdt(n)), 1.0e-30_rt);

            dt = amrex::min(dt, target_dX / dXdt_abs);

        }

    }

    if (integrate_energy) {

        Real target_de;
        if (dedt > 0.0) {
            target_de = (maximum_timestep_change_factor - 1.0_rt) * state.e;
        } else {
            target_de = (1.0_rt - 1.0_rt / maximum_timestep_change_factor) * state.e;
        }

        Real dedt = amrex::max(std::abs(dedt), 1.0e-30_rt);

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

    eos_t eos_state;

    burn_to_eos(state, eos_state);

    eos(eos_input_re, eos_state);

    eos_to_burn(eos_state, state);

    // Ensure that the temperature always stays within reasonable limits.

    state.T = amrex::min(MAX_TEMP, amrex::max(state.T, EOSData::mintemp));

}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void initialize_state (burn_t& state)
{
    // We assume that (rho, T) coming in are valid, do an EOS call
    // to fill the rest of the thermodynamic variables.

    eos_t eos_state;

    burn_to_eos(state, eos_state);

    eos(eos_input_rt, eos_state);

    eos_to_burn(eos_state, state);

    state.self_heat = true;

    state.success = true;

    state.n_rhs = 0;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void evaluate_rhs (burn_t& state, Array1D<Real, 1, NumSpec>& dXdt, Real& dedt)
{
    rate_t rr_old;
    evaluate_rates(state, rr_old);

    rate_fr_t rr;

    // Convert old-style rates to new-style
    for (int j = 1; j <= Rates::NumRatesFR; ++j) {
        rr.forward_rates(j) = rr_old.rates(2 * j - 1);
        rr.reverse_rates(j) = rr_old.rates(2 * j);
    }

    Array1D<Real, 1, NumSpec> y;

    for (int spec = 1; spec <= NumSpec; ++spec) {
        y(spec) = state.xn[spec-1] * aion_inv[spec-1];
    }

    for (int spec = 1; spec <= NumSpec; ++spec) {

        dXdt(spec) = 0.0;

        // Loop through rates.

        for (int rate = 1; rate <= Rates::NumRatesFR; ++rate) {
            rhs_t rhs_data = RHS::rhs_data(rate);

            // First, compute the Y * rate component of both the forward and
            // reverse reactions, which is the same regardless of which species
            // we're producing or consuming.

            Real forward_term = rr.forward_rates(rate);

            if (rhs_data.species_A >= 0) {
                forward_term *= std::pow(y(rhs_data.species_A), rhs_data.exponent_A_forward);
            }

            if (rhs_data.species_B >= 0) {
                forward_term *= std::pow(y(rhs_data.species_B), rhs_data.exponent_B_forward);
            }

            if (rhs_data.species_C >= 0) {
                forward_term *= std::pow(y(rhs_data.species_C), rhs_data.exponent_C_forward);
            }

            Real reverse_term = rr.reverse_rates(rate);

            if (rhs_data.species_D >= 0) {
                reverse_term *= std::pow(y(rhs_data.species_D), rhs_data.exponent_D_reverse);
            }

            if (rhs_data.species_E >= 0) {
                reverse_term *= std::pow(y(rhs_data.species_E), rhs_data.exponent_E_reverse);
            }

            if (rhs_data.species_F >= 0) {
                reverse_term *= std::pow(y(rhs_data.species_F), rhs_data.exponent_F_reverse);
            }

            // Now compute the total contribution to this species.

            if (rhs_data.species_A == spec) {
                forward_term *= -rhs_data.consumed_A_forward;
                reverse_term *= rhs_data.produced_A_reverse;
            }
            else if (rhs_data.species_B == spec) {
                forward_term *= -rhs_data.consumed_B_forward;
                reverse_term *= rhs_data.produced_B_reverse;
            }
            else if (rhs_data.species_C == spec) {
                forward_term *= -rhs_data.consumed_C_forward;
                reverse_term *= rhs_data.produced_C_reverse;
            }
            else if (rhs_data.species_D == spec) {
                forward_term *= rhs_data.produced_D_forward;
                reverse_term *= -rhs_data.consumed_D_reverse;
            }
            else if (rhs_data.species_E == spec) {
                forward_term *= rhs_data.produced_E_forward;
                reverse_term *= -rhs_data.consumed_E_reverse;
            }
            else if (rhs_data.species_F == spec) {
                forward_term *= rhs_data.produced_F_forward;
                reverse_term *= -rhs_data.consumed_F_reverse;
            }
            else {
                forward_term = 0.0_rt;
                reverse_term = 0.0_rt;
            }

            // Scale by aion to get dX/dt instead of dY/dt.
            forward_term *= aion[spec - 1];
            reverse_term *= aion[spec - 1];

            dXdt(spec) += forward_term + reverse_term;
        }
    }

    // Determine the corresponding energy and temperature RHS.

    Array1D<Real, 1, NumSpec> dYdt;
    for (int n = 1; n <= NumSpec; ++n) {
        dYdt(n) = dXdt(n) * aion_inv[n - 1];
    }

    dedt = ener_rhs(state, dYdt);

    state.n_rhs += 1;
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

        Array1D<Real, 1, NumSpec> dXdt;
        Real dedt;

        evaluate_rhs(state, dXdt, dedt);

        state.n_rhs += 1;

        // Calculate the timestep.

        Real dt_sub = calculate_dt(state, dXdt, dedt);

        // Prevent the timestep from overshooting the final time.

        if (t + dt_sub > dt) {
            dt_sub = dt - t;
        }

        for (int n = 1; n <= NumSpec; ++n) {
            state.xn[n-1] += dXdt(n) * dt_sub;
        }

        if (integrate_energy) {
            state.e += dedt * dt_sub;
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
        std::cout <<  "energy released: " << state.e - e_in << std::endl;
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
