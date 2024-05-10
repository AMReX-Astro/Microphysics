#ifndef actual_integrator_H
#define actual_integrator_H

// Common variables and routines for burners
// that use VODE for their integration.

#include <AMReX_Print.H>

#include <iomanip>

#include <network.H>
#include <burn_type.H>
#include <eos_type.H>
#include <eos.H>
#include <extern_parameters.H>
#include <vode_type.H>
#include <vode_dvode.H>

using namespace integrator_rp;

template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_integrator (BurnT& state, amrex::Real dt, bool is_retry=false)
{
    constexpr int int_neqs = integrator_neqs<BurnT>();

    dvode_t<int_neqs> vode_state{};

    // Start off by assuming a successful burn.

    state.success = true;

    // Initialize the integration time.

    vode_state.t = 0.0_rt;
    vode_state.tout = dt;

    // Set the (inverse of the) timestep limiter.

    vode_state.HMXI = 1.0_rt / ode_max_dt;

    // set the Jacobian type
    if (is_retry && integrator_rp::retry_swap_jacobian) {
        vode_state.jacobian_type = (jacobian == 1) ? 2 : 1;
    } else {
        vode_state.jacobian_type = jacobian;
    }

    // Fill in the initial integration state.

    burn_to_int(state, vode_state);

    // Save the initial composition and temperature for our later diagnostics.

#ifndef AMREX_USE_GPU
    amrex::Real xn_in[NumSpec];
    for (int n = 0; n < NumSpec; ++n) {
        xn_in[n] = state.y[SFS+n] / state.y[SRHO];
    }
    // we are assuming that the temperature was valid on input
    amrex::Real T_in = state.T;
#ifdef AUX_THERMO
    amrex::Real aux_in[NumAux];
    for (int n = 0; n < NumAux; ++n) {
        aux_in[n] = state.y[SFX+n] / state.y[SRHO];
    }
#endif
    amrex::Real rhoe_in = state.y[SEINT];
#endif


    // Set the tolerances.

    amrex::Real sdc_tol_fac = std::pow(integrator_rp::sdc_burn_tol_factor, state.num_sdc_iters - state.sdc_iter - 1);

    // we use 1-based indexing inside of VODE, so we need to shift the
    // indices SRHO, SFS, etc by 1

    amrex::Real sdc_min_density = amrex::min(state.rho, state.rho_orig + state.ydot_a[SRHO] * dt);

    if (!is_retry) {

        vode_state.atol_enuc = sdc_min_density * atol_enuc * sdc_tol_fac;
        vode_state.rtol_enuc = rtol_enuc * sdc_tol_fac;

        // Note: we define the input atol for species to refer only to the
        // mass fraction part, and we multiply by a representative density
        // so that atol becomes an absolutely tolerance on (rho X)

        vode_state.atol_spec = sdc_min_density * atol_spec * sdc_tol_fac;
        vode_state.rtol_spec = rtol_spec * sdc_tol_fac;

    } else {

        vode_state.atol_enuc = sdc_min_density * retry_atol_enuc * sdc_tol_fac;
        vode_state.rtol_enuc = retry_rtol_enuc * sdc_tol_fac;

        // Note: we define the input atol for species to refer only to the
        // mass fraction part, and we multiply by a representative density
        // so that atol becomes an absolutely tolerance on (rho X)

        vode_state.atol_spec = sdc_min_density * retry_atol_spec * sdc_tol_fac;
        vode_state.rtol_spec = retry_rtol_spec * sdc_tol_fac;

    }

    if (scale_system) {
        // the absolute tol for energy needs to reflect the scaled
        // energy the integrator sees
        vode_state.atol_enuc /= state.e_scale;
    }

    // Call the integration routine.

    int istate = dvode(state, vode_state);
    state.error_code = istate;

    // Get the number of RHS and Jacobian evaluations.

    state.n_rhs = vode_state.n_rhs;
    state.n_jac = vode_state.n_jac;
    state.n_step = vode_state.n_step;

    // Copy the integration data back to the burn state.
    // This will also update the aux state from X if we are using NSE

    int_to_burn(vode_state.t, vode_state, state);

    // we only evolved (rho e), not (rho E), so we need to update the
    // total energy now to ensure we are conservative

    amrex::Real rho_Sdot = 0.0_rt;
    if (state.time > 0) {
        rho_Sdot = (state.y[SEINT] - state.rhoe_orig) / state.time - state.ydot_a[SEINT];
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

    // VODE does not always fail even though it can lead to unphysical states.
    // Add some checks that indicate a burn fail even if VODE thinks the
    // integration was successful.

    if (istate != IERR_SUCCESS) {
        state.success = false;
    }

    if (state.y[SEINT] < 0.0_rt) {
        state.success = false;
    }

    for (int n = 0; n < NumSpec; ++n) {
        if (state.y[SFS+n] / state.rho < -species_failure_tolerance) {
            state.success = false;
        }

        if (state.y[SFS+n] / state.rho > 1.0_rt + species_failure_tolerance) {
            state.success = false;
        }
    }


#ifndef AMREX_USE_GPU
    if (burner_verbose) {
        // Print out some integration statistics, if desired.
        std::cout <<  "integration summary: " << std::endl;
        std::cout <<  "dens: " << state.rho << " temp: " << state.T << std::endl;
        std::cout << " energy released: " << state.e << std::endl;
        std::cout <<  "number of steps taken: " << state.n_step << std::endl;
        std::cout <<  "number of f evaluations: " << state.n_rhs << std::endl;
    }
#endif

    // If we failed, print out the current state of the integration.

    if (!state.success) {
        if (istate != IERR_ENTERED_NSE) {
#ifndef AMREX_USE_GPU
            std::cout << amrex::Font::Bold << amrex::FGColor::Red << "[ERROR] integration failed in net" << amrex::ResetDisplay << std::endl;
            std::cout << "istate = " << istate << std::endl;
            if (istate == IERR_SUCCESS) {
                std::cout << "  VODE exited successfully, but a check on the data values failed" << std::endl;
            }
            std::cout << "zone = (" << state.i << ", " << state.j << ", " << state.k << ")" << std::endl;
            std::cout << "time = " << state.time << std::endl;
            std::cout << "dt = " << std::setprecision(16) << dt << std::endl;
            std::cout << "dens start = " << std::setprecision(16) << state.rho_orig << std::endl;
            std::cout << "temp start = " << std::setprecision(16) << T_in << std::endl;
            std::cout << "rhoe start = " << std::setprecision(16) << rhoe_in << std::endl;
            std::cout << "xn start = ";
            for (const auto X : xn_in) {
                std::cout << std::setprecision(16) << X << " ";
            }
            std::cout << std::endl;
#ifdef AUX_THERMO
            std::cout << "aux start = ";
            for (const auto aux : aux_in) {
                std::cout << std::setprecision(16) << aux << " ";
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
        } else {
#ifndef AMREX_USE_GPU
            std::cout << "burn entered NSE during integration (after " << state.n_step << " steps), zone = (" << state.i << ", " << state.j << ", " << state.k << ")" << std::endl;
#endif
        }
    }

}

#endif
