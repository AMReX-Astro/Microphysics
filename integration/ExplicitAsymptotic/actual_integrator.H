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
void clean_state (burn_t& state);

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void initialize_state (burn_t& state);

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void evaluate_rhs (burn_t& state, Array1D<Real, 1, NumSpec>& f_minus, Array1D<Real, 1, NumSpec>& f_plus,
                   Real& dedt, Real& dTdt);

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void calculate_dydt (burn_t& state, const Real& dt,
                     const Array1D<Real, 1, NumSpec>& f_minus, const Array1D<Real, 1, NumSpec>& f_plus,
                     Array1D<Real, 1, NumSpec>& dXdt, Real& dedt, Real& dTdt);

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real limit_dt (const burn_t& state, const Array1D<Real, 1, NumSpec>& dXdt, const Real& dt);

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void actual_integrator (burn_t& state, Real dt)
{
    initialize_state(state);

    Real t = 0.0;

    // Start the guess for the timestepping with a specified fraction of the full dt.

    Real dt_sub = dt * dt_init_fraction;

    // When checking the integration time to see if we're done,
    // be careful with roundoff issues.

    const Real timestep_safety_factor = 1.0e-12_rt;

    int num_timesteps = 0;

    while (t < (1.0_rt - timestep_safety_factor) * dt)
    {
        // Evaluate the RHS.

        Array1D<Real, 1, NumSpec> f_plus, f_minus;
        Real dedt, dTdt;

        evaluate_rhs(state, f_minus, f_plus, dedt, dTdt);

        // Iterate over dt, first calculating the update and then determining
        // whether the update satisfies any constraints we have on the timestep.
        // Our initial guess is the previous timestep multiplied by a factor > 1
        // to allow the timestep to grow in the absence of any constraints.

        Real dt_new = dt_max_change_factor * dt_sub;
        int num_iters = 0;

        Array1D<Real, 1, NumSpec> dXdt;

        while (std::abs(dt_new - dt_sub) > dt_error_tolerance * dt_sub && num_iters < dt_max_iters)
        {
            dt_sub = dt_new;
            ++num_iters;

            calculate_dydt(state, dt_sub, f_minus, f_plus, dXdt, dedt, dTdt);

            dt_new = limit_dt(state, dXdt, dt_sub);
        }

        // Update the timestep from the last iteration.

        dt_sub = dt_new;

        // Prevent the timestep from overshooting the final time.

        if (t + dt_sub > dt) {
            dt_sub = dt - t;
        }

        // Update the state quantities.

        for (int n = 1; n <= NumSpec; ++n) {
            state.xn[n-1] += dXdt(n) * dt_sub;
        }

        if (integrate_temperature) {
            state.T += dTdt * dt_sub;
        }
        if (integrate_energy) {
            state.e += dedt * dt_sub;
        }

        clean_state(state);

        t += dt_sub;
        ++num_timesteps;
    }

#ifndef AMREX_USE_CUDA
    if (burner_verbose) {
        // Print out some integration statistics, if desired.
        std::cout <<  "integration summary: " << std::endl;
        std::cout <<  "dens: " << state.rho << " temp: " << state.T << std::endl;
        std::cout <<  "energy released: " << state.e << std::endl;
        std::cout <<  "number of steps taken: " << num_timesteps << std::endl;
        std::cout <<  "number of f evaluations: " << state.n_rhs << std::endl;
    }
#endif
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
void initialize_state (burn_t& state)
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
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void evaluate_rhs (burn_t& state, Array1D<Real, 1, NumSpec>& f_minus, Array1D<Real, 1, NumSpec>& f_plus,
                   Real& dedt, Real& dTdt)
{
    rate_t rr;
    evaluate_rates(state, rr);

    Array1D<Real, 1, NumSpec> y;

    for (int spec = 1; spec <= NumSpec; ++spec) {
        y(spec) = state.xn[spec-1] * aion_inv[spec-1];
    }

    // For each species calculate its positive and negative RHS contributions.

    for (int spec = 1; spec <= NumSpec; ++spec) {

        f_plus(spec) = 0.0;
        f_minus(spec) = 0.0;

        // Loop through rates. All terms that have a positive contribution
        // to the RHS for this species get the full term added to f_plus.
        // All terms that deplete this species are added to f_minus.

        for (int rate = 1; rate <= Rates::NumRates; ++rate) {
            rhs_t rhs_data = RHS::rhs_data(spec, rate);

            if (std::abs(rhs_data.prefactor) > 0.0) {

                Real term = rhs_data.prefactor;

                if (rhs_data.specindex1 >= 0) {
                    term *= y(rhs_data.specindex1);
                }

                if (rhs_data.specindex2 >= 0) {
                    term *= y(rhs_data.specindex2);
                }

                if (rhs_data.specindex3 >= 0) {
                    term *= y(rhs_data.specindex3);
                }

                term *= rr.rates(rate);

                // Scale by aion to get dX/dt instead of dY/dt.

                term *= aion[spec - 1];

                if (rhs_data.prefactor > 0.0) {
                    f_plus(spec) += term;
                } else {
                    f_minus(spec) -= term;
                }
            }
        }
    }
 
    // Determine the corresponding energy and temperature RHS.

    Array1D<Real, 1, NumSpec> dYdt;
    for (int n = 1; n <= NumSpec; ++n) {
        dYdt(n) = (f_plus(n) - f_minus(n)) * aion_inv[n - 1];
    }

    dedt = ener_rhs(state, dYdt);
    dTdt = temperature_rhs(state, dedt);

    state.n_rhs += 1;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void calculate_dydt (burn_t& state, const Real& dt,
                     const Array1D<Real, 1, NumSpec>& f_minus, const Array1D<Real, 1, NumSpec>& f_plus,
                     Array1D<Real, 1, NumSpec>& dXdt, Real& dedt, Real& dTdt)
{
    // Given our guess for dt, determine whether this species should be updated using
    // the asymptotic approximation. The decision criterion for allowing the
    // asymptotic approximation is whether (f_minus / X) * dt >= 1 (section 7.3).
    // If we're using the approximation, we calculate the new-time value of X
    // using Equation 13, and then calculate an effective "explicit" dX/dt based on
    // this, so that we can just update everything by dt * dX/dt at the end.

    for (int n = 1; n < NumSpec; ++n)
    {
        if (f_minus(n) * dt >= state.xn[n-1])
        {
            Real X_new = (state.xn[n-1] + f_plus(n) * dt) / (1 + f_minus(n) * dt / state.xn[n-1]);
            dXdt(n) = (X_new - state.xn[n-1]) / dt;
        }
        else
        {
            dXdt(n) = (f_plus(n) - f_minus(n));
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real limit_dt (const burn_t& state, const Array1D<Real, 1, NumSpec>& dXdt, const Real& dt)
{
    // Limit dt by whether the update would change the abundance
    // by more than a given factor, using simple linear limiting.
    // This step is skipped if the species is below atol_spec.

    Real dt_new = dt;

    for (int n = 1; n < NumSpec; ++n) {
        if (state.xn[n-1] >= atol_spec) {
            dt_new = amrex::min(dt_new, std::abs((maximum_X_change_factor_per_timestep - 1.0_rt) * state.xn[n-1] / dXdt(n)));
        }
    }

    dt_new = amrex::min(dt_new, ode_max_dt);

    return dt_new;
}

#endif
