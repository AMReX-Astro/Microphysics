#ifndef actual_integrator_H
#define actual_integrator_H

#include <iomanip>

#include <be_type.H>
#include <be_integrator.H>
#include <extern_parameters.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void actual_integrator (burn_t& state, const Real dt)
{

    be_t be;

    // set the Jacobian type
    be.jacobian_type = jacobian;

    // Start off by assuming a successful burn.

    state.success = true;

    // Initialize the integration time.

    be.t = 0.0;
    be.tout = dt;

    // Fill in the initial integration state.

    burn_to_int(state, be);

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
    Real rhoe_in = state.y[SEINT];
#endif


    // Set the tolerances.

    Real sdc_tol_fac = std::pow(sdc_burn_tol_factor, state.num_sdc_iters - state.sdc_iter - 1);

    // we use 1-based indexing inside of BackwardEuler, so we need to shift the
    // indices SRHO, SFS, etc by 1

    Real sdc_min_density = amrex::min(state.rho, state.rho_orig + state.ydot_a[SRHO] * dt);

    be.atol_enuc = sdc_min_density * atol_enuc * sdc_tol_fac;
    be.rtol_enuc = rtol_enuc * sdc_tol_fac;

    // Note: we define the input atol for species to refer only to the
    // mass fraction part, and we multiply by a representative density
    // so that atol becomes an absolutely tolerance on (rho X)

    be.atol_spec = sdc_min_density * atol_spec * sdc_tol_fac;
    be.rtol_spec = rtol_spec * sdc_tol_fac;

    // Call the integration routine.

    int istate = be_integrator(state, be);

    // Get the number of RHS and Jacobian evaluations.

    state.n_rhs = be.n_rhs;
    state.n_jac = be.n_jac;
    state.n_step = be.n_step;

    // Copy the integration data back to the burn state.
    // This will also update the aux state from X if we are using NSE

    int_to_burn(be.t, be, state);

    // we only evolved (rho e), not (rho E), so we need to update the
    // total energy now to ensure we are conservative

    Real rho_Sdot = 0.0_rt;
    if (state.time > 0) {
        rho_Sdot = (be.y(SEINT+1) - state.rhoe_orig) / state.time - state.ydot_a[SEINT];
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

    if (istate < 0) {
        state.success = false;
    }


#ifndef AMREX_USE_GPU
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
#ifndef AMREX_USE_GPU
        std::cout << Font::Bold << FGColor::Red << "[ERROR] integration failed in net" << ResetDisplay << std::endl;
        std::cout << "istate = " << istate << std::endl;
        std::cout << "zone = (" << state.i << ", " << state.j << ", " << state.k << ")" << std::endl;
        std::cout << "time = " << be.t << std::endl;
        std::cout << "dt = " << std::setprecision(16) << dt << std::endl;
        std::cout << "dens start = " << std::setprecision(16) << state.rho_orig << std::endl;
        std::cout << "temp start = " << std::setprecision(16) << T_in << std::endl;
        std::cout << "rhoe start = " << std::setprecision(16) << rhoe_in << std::endl;
        std::cout << "xn start = ";
        for (int n = 0; n < NumSpec; ++n) {
            std::cout << std::setprecision(16) << xn_in[n] << " ";
        }
        std::cout << std::endl;
#ifdef AUX_THERMO
        std::cout << "aux start = ";
        for (int n = 0; n < NumAux; ++n) {
            std::cout << std::setprecision(16) << aux_in[n] << " ";
        }
        std::cout << std::endl;
#endif
        std::cout << "dens current = " << std::setprecision(16) << state.rho << std::endl;
        std::cout << "temp current = " << std::setprecision(16) << state.T << std::endl;
        std::cout << "xn current = ";
        for (int n = 0; n < NumSpec; ++n) {
            std::cout << std::setprecision(16) << state.xn[n] << " ";
        }
        std::cout << std::endl;
#ifdef AUX_THERMO
        std::cout << "aux current = ";
        for (int n = 0; n < NumAux; ++n) {
            std::cout << std::setprecision(16) << state.aux[n] << " ";
        }
        std::cout << std::endl;
#endif
        std::cout << "A(rho) = " << std::setprecision(16) << state.ydot_a[SRHO] << std::endl;
        std::cout << "A(rho e) = " << std::setprecision(16) << state.ydot_a[SEINT] << std::endl;
        std::cout << "A(rho X_k) = ";
        for (int n = 0; n < NumSpec; n++) {
            std::cout << std::setprecision(16) << state.ydot_a[SFS+n] << " ";
        }
        std::cout << std::endl;
#ifdef AUX_THERMO
        std::cout << "A(rho aux_k) = ";
        for (int n = 0; n < NumAux; n++) {
            std::cout << std::setprecision(16) << state.ydot_a[SFX+n] << " ";
        }
        std::cout << std::endl;
#endif
#endif
    }


}

#endif
