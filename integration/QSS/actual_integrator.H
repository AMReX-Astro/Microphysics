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
                   Real& dedt);

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void calculate_dydt (burn_t& state, const Real& dt, Array1D<Real, 1, NumSpec>& dXdt, Real& dedt);

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real limit_dt (const burn_t& state, const Array1D<Real, 1, NumSpec>& dXdt, const Real& dedt, const Real& dt);

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

    // Start the guess for the timestepping with a specified fraction of the full dt.

    Real dt_sub = dt * dt_init_fraction;

    // When checking the integration time to see if we're done,
    // be careful with roundoff issues.

    const Real timestep_safety_factor = 1.0e-12_rt;

    int num_timesteps = 0;

    while (t < (1.0_rt - timestep_safety_factor) * dt && num_timesteps < ode_max_steps)
    {
        // Iterate over dt, first calculating the update and then determining
        // whether the update satisfies any constraints we have on the timestep.
        // Our initial guess is the previous timestep multiplied by a factor > 1
        // to allow the timestep to grow in the absence of any constraints.

        Real dt_new = dt_max_change_factor * dt_sub;
        int num_iters = 0;

        Array1D<Real, 1, NumSpec> dXdt;
        Real dedt;

        while (std::abs(dt_new - dt_sub) > dt_error_tolerance * dt_sub && num_iters < dt_max_iters)
        {
            dt_sub = dt_new;
            ++num_iters;

            calculate_dydt(state, dt_sub, dXdt, dedt);

            dt_new = limit_dt(state, dXdt, dedt, dt_sub);
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

        state.e += dedt * dt_sub;

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
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void evaluate_rhs (burn_t& state, Array1D<Real, 1, NumSpec>& f_minus, Array1D<Real, 1, NumSpec>& f_plus,
                   Real& dedt)
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

    state.n_rhs += 1;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real alpha (const Real& r)
{
    // Evaluate alpha (Guidry et al. paper II, Equation 4)
    // Note the first term in the numerator has a typo,
    // compare to Equation 39 of Mott et al. 2000.

    if (std::abs(r) <= 1.e-30_rt) { // Protect against underflow
        return 1.0_rt; // Limit for r -> 0
    }
    else if (std::abs(r) >= 1.e30_rt) { // Protect against overflow
        return 0.5_rt; // Limit for r -> infinity
    }
    else {
        Real r2 = r * r;
        Real r3 = r2 * r;
        Real a = (180.0_rt * r3 + 60.0_rt * r2 + 11.0_rt * r + 1.0_rt);
        a /= (360.0_rt * r3 + 60.0_rt * r2 + 12.0_rt * r + 1.0_rt);
        return a;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void calculate_dydt (burn_t& state, const Real& dt, Array1D<Real, 1, NumSpec>& dXdt, Real& dedt)
{
    // Calculate the quasi-steady-state update. This is based on Guidry et al. paper II.
    // Section 2 spells out the update for the species using a predictor-corrector method.
    // The RHS for any term is split into terms f_plus and f_minus, where generally the
    // terms in f_plus add to X, while f_minus contains all terms that deplete X.
    // (In practice, we have some two-body reactions like C12 + C12 and three-body reactions
    // like triple alpha where k is a function of X. The method ignores this subtlety.)
    //
    // dX/dt = f_plus - f_minus = f_plus - k * X
    //
    // In the QSS method we define the following variables:
    //
    // r = 1 / (k * dt)
    // alpha(r) = (180 * r**3 + 60 * r**2 + 11 * r + 1) / (360 * r**3 + 60 * r**2 + 12 * r + 1)
    //
    // and then do a predictor step defined by Guidry Equation 5
    //
    // X_p = X_0 + dt * (f_plus - k_0 * X_0) / (1 + alpha_0 * k_0 * dt)
    //
    // Then we recompute f_plus and f_minus (and therefore k) from the RHS using the predictor state,
    // and compute time-centered variables
    //
    // k_bar = (k_0 + k_p) / 2
    // alpha_bar = (alpha_0 + alpha_p) / 2
    //
    // and then recompute the corrector state
    //
    // X_c = X_0 + (f_plus_tilde - k_bar * X_0) / (1 + alpha_bar * k_bar * dt)
    //
    // where f_plus_tilde is a weighted average of the old and new terms:
    //
    // f_plus_tilde = alpha_bar * f_plus_p + (1 - alpha) * f_plus_0
    //
    // We can also do predictor-corrector on the energy update, which has the more familiar form
    //
    // e_p = e_0 + dt * f_0
    // e_c = e_0 + dt * (f_0 + f_p) / 2
    //
    // This is actually a special case of the QSS update for k == 0 (no "self-depletion" term) and
    // alpha = 1/2.
    //
    // We do not integrate temperature, and instead get it from the EOS.

    // Make a copy of the old-time state; we'll be temporarily updating the state with the predictor
    // values for the new-time RHS evaluation, but we need to revert those changes later.

    burn_t state_0 = state;

    // Evaluate the RHS using the old-time state.

    Array1D<Real, 1, NumSpec> f_plus_0, f_minus_0;
    Real dedt_0;

    evaluate_rhs(state_0, f_minus_0, f_plus_0, dedt_0);

    for (int n = 1; n <= NumSpec; ++n)
    {
        Real X_0 = state_0.xn[n-1];
        Real k_0 = f_minus_0(n) / X_0;
        Real r_0 = 1.0_rt / amrex::max(k_0 * dt, 1.0e-50_rt);
        Real alpha_0 = alpha(r_0);

        state.xn[n-1] = X_0 + dt * (f_plus_0(n) - k_0 * X_0) / (1.0_rt + alpha_0 * k_0 * dt);
    }

    state.e = state_0.e + dt * dedt_0;

    // Now that we have the predictor, re-evaluate the RHS. Note that this also
    // involves a re-evaluation of the rates. This seems safe given how temperature
    // sensitive our rates are. Mott et al. recommend this for cases where the integration
    // is sensitive to T (second to last paragraph of section 4).

    clean_state(state);

    Array1D<Real, 1, NumSpec> f_minus_p, f_plus_p;
    Real dedt_p;

    evaluate_rhs(state, f_minus_p, f_plus_p, dedt_p);

    // Compute the corrector state.

    for (int n = 1; n <= NumSpec; ++n)
    {
        Real X_0 = state_0.xn[n-1];
        Real X_p = state.xn[n-1];
        Real k_0 = f_minus_0(n) / X_0;
        Real k_p = f_minus_p(n) / X_p;
        Real k_bar = 0.5_rt * (k_p + k_0);
        Real r_bar = 1.0_rt / amrex::max(k_bar * dt, 1.0e-50_rt);
        Real alpha_bar = alpha(r_bar);
        Real f_plus_tilde = alpha_bar * f_plus_p(n) + (1.0_rt - alpha_bar) * f_plus_0(n);

        state.xn[n-1] = X_0 + dt * (f_plus_tilde - k_bar * X_0) / (1.0_rt + alpha_bar * k_bar * dt);
    }

    state.e = state_0.e + 0.5_rt * dt * (dedt_0 + dedt_p);

    // Now compute the final values for dX/dt and de/dt, and then
    // reset the state back to the way it was before.

    for (int n = 1; n <= NumSpec; ++n)
    {
        dXdt(n) = (state.xn[n-1] - state_0.xn[n-1]) / dt;
    }

    dedt = (state.e - state_0.e) / dt;

    state = state_0;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real limit_dt (const burn_t& state, const Array1D<Real, 1, NumSpec>& dXdt, const Real& dedt, const Real& dt)
{
    // Limit dt by whether the update would change the abundance
    // by more than a given factor, using simple linear limiting.
    // This step is skipped if the species is below atol_spec.

    Real dt_new = dt;

    for (int n = 1; n < NumSpec; ++n) {
        if (state.xn[n-1] >= atol_spec) {
            dt_new = amrex::min(dt_new, std::abs((maximum_change_factor_per_timestep - 1.0_rt) * state.xn[n-1] / dXdt(n)));
        }
    }

    dt_new = amrex::min(dt_new, std::abs((maximum_change_factor_per_timestep - 1.0_rt) * state.e / dedt));

    dt_new = amrex::min(dt_new, ode_max_dt);

    return dt_new;
}

#endif
