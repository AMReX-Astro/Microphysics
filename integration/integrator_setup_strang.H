#ifndef INTEGRATOR_SETUP_STRANG_H
#define INTEGRATOR_SETUP_STRANG_H

#include <AMReX_REAL.H>
#include <AMReX_ANSIEscCode.H>

#include <iomanip>

#include <network.H>
#include <burn_type.H>
#include <eos_type.H>
#include <eos.H>
#include <integrator_data.H>

#include <extern_parameters.H>

struct state_backup_t {
    amrex::Real T_in{};
    amrex::Real e_in{};
#ifndef AMREX_USE_GPU
    amrex::Real xn_in[NumSpec]{};
#endif
};

template <typename BurnT, typename IntegratorT>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
IntegratorT integrator_setup (BurnT& state, amrex::Real dt, bool is_retry)
{

    IntegratorT int_state{};

    // Set the tolerances.

    if (!is_retry) {
        int_state.atol_spec = integrator_rp::atol_spec; // mass fractions
        int_state.atol_enuc = integrator_rp::atol_enuc; // energy generated

        int_state.rtol_spec = integrator_rp::rtol_spec; // mass fractions
        int_state.rtol_enuc = integrator_rp::rtol_enuc; // energy generated
    } else {
        int_state.atol_spec = integrator_rp::retry_atol_spec; // mass fractions
        int_state.atol_enuc = integrator_rp::retry_atol_enuc; // energy generated

        int_state.rtol_spec = integrator_rp::retry_rtol_spec; // mass fractions
        int_state.rtol_enuc = integrator_rp::retry_rtol_enuc; // energy generated
    }

    // set the Jacobian type
    if (is_retry && integrator_rp::retry_swap_jacobian) {
        int_state.jacobian_type = (integrator_rp::jacobian == 1) ? 2 : 1;
    } else {
        int_state.jacobian_type = static_cast<short>(integrator_rp::jacobian);
    }

    // Start off by assuming a successful burn.

    state.success = true;

    // Initialize the integration time.

    int_state.t = 0.0_rt;
    int_state.tout = dt;

    // Initialize ydot to zero for Strang burn.

    for (int n = 0; n < SVAR; ++n) {
        state.ydot_a[n] = 0;
    }

    // We assume that (rho, T) coming in are valid, do an EOS call
    // to fill the rest of the thermodynamic variables.

    eos(eos_input_rt, state);

    // set the scaling for energy if we integrate it dimensionlessly
    state.e_scale = state.e;

    if (integrator_rp::scale_system) {
        // the absolute tol for energy needs to reflect the scaled
        // energy the integrator sees
        int_state.atol_enuc /= state.e_scale;
    }

    // Fill in the initial integration state.

    burn_to_integrator(state, int_state);

    // Save the initial composition, temperature, and energy for our later diagnostics.

    return int_state;
}


template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
state_backup_t integrator_backup (const BurnT& state) {

    state_backup_t state_save;

#ifndef AMREX_USE_GPU
    for (int n = 0; n < NumSpec; ++n) {
        state_save.xn_in[n] = state.xn[n];
    }
    state_save.T_in = state.T;
#endif
    state_save.e_in = state.e;

    return state_save;

}


template <typename BurnT, typename IntegratorT>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void integrator_cleanup (IntegratorT& int_state, BurnT& state,
                         int istate, const state_backup_t& state_save, amrex::Real dt)
{

    // Copy the integration data back to the burn state.

    integrator_to_burn(int_state, state);

#ifdef NSE
    // compute the temperature based on the energy release -- we need
    // this in case we failed in our burn here because we entered NSE

#ifdef AUX_THERMO
    // need to sync the auxiliary data up with the new mass fractions
    set_aux_comp_from_X(state);
#endif
    if (integrator_rp::call_eos_in_rhs) {
        eos(eos_input_re, state);
    }
#endif

    // Subtract off the initial energy (the application codes expect
    // to get back only the generated energy during the burn).

    if (integrator_rp::subtract_internal_energy) {
        state.e -= state_save.e_in;
    }

    // Normalize the final abundances (except if they are number
    // densities)

    if (! integrator_rp::use_number_densities) {
        normalize_abundances_burn(state);
    }

    // Get the number of RHS and Jacobian evaluations.

    state.n_rhs = int_state.n_rhs;
    state.n_jac = int_state.n_jac;
    state.n_step = int_state.n_step;

    // The integrator may not always fail even though it can lead to
    // unphysical states.  Add some checks that indicate a burn fail
    // even if the integrator thinks the integration was successful.

    if (istate != IERR_SUCCESS) {
        state.success = false;
    }

    for (int n = 1; n <= NumSpec; ++n) {
        if (int_state.y(n) < -species_failure_tolerance) {
            state.success = false;
        }

        // Don't enforce a max if we are evolving number densities

        if (! integrator_rp::use_number_densities) {
            if (int_state.y(n) > 1.0_rt + species_failure_tolerance) {
                state.success = false;
            }
        }
    }

#ifndef AMREX_USE_GPU
    if (integrator_rp::burner_verbose) {
        // Print out some integration statistics, if desired.
        std::cout <<  "integration summary: " << std::endl;
        std::cout <<  "dens: " << state.rho << " temp: " << state.T << std::endl;
        std::cout << " energy released: " << state.e << std::endl;
        std::cout <<  "number of steps taken: " << state.n_step << std::endl;
        std::cout <<  "number of f evaluations: " << state.n_rhs << std::endl;
    }
#endif

    // If we failed, print out the current state of the integration.

    if (! state.success) {
        if (istate != IERR_ENTERED_NSE) {
#ifndef AMREX_USE_GPU
            std::cout << amrex::Font::Bold << amrex::FGColor::Red << "[ERROR] integration failed in net" << amrex::ResetDisplay << std::endl;
            std::cout << "istate = " << istate << std::endl;
            if (istate == IERR_SUCCESS) {
                std::cout << "  Integration exited successfully, but a check on the data values failed" << std::endl;
            }
            std::cout << "zone = (" << state.i << ", " << state.j << ", " << state.k << ")" << std::endl;
            std::cout << "time = " << int_state.t << std::endl;
            std::cout << "dt = " << std::setprecision(16) << dt << std::endl;
            std::cout << "temp start = " << std::setprecision(16) << state_save.T_in << std::endl;
            std::cout << "xn start = ";
            for (const double X : state_save.xn_in) {
                std::cout << std::setprecision(16) << X << " ";
            }
            std::cout << std::endl;
            std::cout << "dens current = " << std::setprecision(16) << state.rho << std::endl;
            std::cout << "temp current = " << std::setprecision(16) << state.T << std::endl;
            std::cout << "xn current = ";
            for (const double X : state.xn) {
                std::cout << std::setprecision(16) << X << " ";
            }
            std::cout << std::endl;
            std::cout << "energy generated = " << state.e << std::endl;
#endif
        } else {
#ifndef AMREX_USE_GPU
            std::cout << "burn entered NSE during integration (after " << state.n_step << " steps), zone = (" << state.i << ", " << state.j << ", " << state.k << ")" << std::endl;
#endif
        }
    }

}

#endif
