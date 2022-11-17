#ifndef actual_integrator_H
#define actual_integrator_H

#include <be_integrator.H>
#include <be_type_strang.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void actual_integrator (burn_t& state, const Real dt)
{

    // start off by assuming a successful burn

    state.success = true;

    be_t be;

    // set the Jacobian type
    be.jacobian_type = jacobian;

    // Initialize the integration time.
    be.t = 0.0;
    be.tout = dt;

    // Fill the initial integration state

    burn_to_be(state, be);

    // Save the initial composition and temperature for our later diagnostics.

#ifndef AMREX_USE_GPU
    Real xn_in[NumSpec];
    for (int n = 0; n < NumSpec; ++n) {
        xn_in[n] = state.y[SFS+n] / state.y[SRHO];
    }
    // we are assuming that the temperature was valid on input
    Real T_in = state.T;
#ifdef AUX_THERMO
    Real aux_in[NumAux];
    for (int n = 0; n < NumAux; ++n) {
        aux_in[n] = state.y[SFX+n] / state.y[SRHO];
    }
#endif
#endif

    // Set the tolerances

    Real sdc_tol_fac = std::pow(sdc_burn_tol_factor, state.num_sdc_iters - state.sdc_iter - 1);

    Real sdc_min_density = amrex::min(state.rho, state.rho + state.ydot_a[SRHO] * dt);

    be.atol_enuc = sdc_min_density * atol_enuc * sdc_tol_fac;
    be.rtol_enuc = rtol_enuc * sdc_tol_fac;

    // Note: we define the input atol for species to refer only to the
    // mass fraction part, and we multiply by a representative density
    // so that atol becomes an absolutely tolerance on (rho X)

    be.atol_spec = sdc_min_density * atol_spec * sdc_tol_fac;
    be.rtol_spec = rtol_spec * sdc_tol_fac;

    // Call the integration routine.

    int ierr = be_integrator(state, be);

    // Get the number of RHS and Jacobian evaluations.

    state.n_rhs = be.n_rhs;
    state.n_jac = be.n_jac;
    state.n_step = be.n_step;

    // Copy the integration data back to the burn state.
    // This will also update the aux state from X if we are using NSE

    be_to_burn(be.t, be, state);

    // we only evolved (rho e), not (rho E), so we need to update the
    // total energy now to ensure we are conservative

    Real rho_Sdot = 0.0_rt;
    if (state.time > 0) {
        rho_Sdot = (vode_state.y(SEINT+1) - state.rhoe_orig) / state.time - state.ydot_a[SEINT];
    }

    state.y[SEDEN] += state.time * (state.ydot_a[SEDEN] + rho_Sdot);

    // also momentum

    state.y[SMX] += state.time * state.ydot_a[SMX];
    state.y[SMY] += state.time * state.ydot_a[SMY];
    state.y[SMZ] += state.time * state.ydot_a[SMZ];

    // normalize the abundances on exit.  We'll assume that the driver
    // calling this is making use of the conserved state (state.y[]),
    // so that is what will be normalized.

    normalize_abundances_sdc_burn(state);

    if (ierr < 0) {
        state.success = false;
    }


#ifndef AMREX_USE_CUDA
    if (burner_verbose) {
        // Print out some integration statistics, if desired.
        std::cout <<  "integration summary: " << std::endl;
        std::cout <<  "dens: " << state.rho << " temp: " << state.T << std::endl;
        std::cout <<  "energy released: " << state.e << std::endl;
        std::cout <<  "number of steps taken: " << state.n_step << std::endl;
        std::cout <<  "number of f evaluations: " << state.n_rhs << std::endl;
    }
#endif

    // If we failed, print out the current state of the integration.

    if (!state.success) {
#ifndef AMREX_USE_CUDA
        std::cout << "ERROR: integration failed in net" << std::endl;
        std::cout << "time = " << be.t << std::endl;
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
    }


}

#endif
