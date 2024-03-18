#ifndef actual_integrator_H
#define actual_integrator_H

#include <iomanip>

#include <network.H>
#include <burn_type.H>

#include <be_type.H>
#include <be_integrator.H>

template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_integrator (BurnT& state, const amrex::Real dt, bool is_retry=false)
{
    constexpr int int_neqs = integrator_neqs<BurnT>();

    be_t<int_neqs> be;

    // Set the tolerances.

    if (!is_retry) {
        be.atol_spec = atol_spec;  // mass fractions
        be.atol_enuc = atol_enuc;  // energy generated

        be.rtol_spec = rtol_spec;  // mass fractions
        be.rtol_enuc = rtol_enuc;  // energy generated
    } else {
        be.atol_spec = retry_atol_spec; // mass fractions
        be.atol_enuc = retry_atol_enuc; // energy generated

        be.rtol_spec = retry_rtol_spec; // mass fractions
        be.rtol_enuc = retry_rtol_enuc; // energy generated
    }

    // set the Jacobian type
    if (is_retry && integrator_rp::retry_swap_jacobian) {
        be.jacobian_type = (jacobian == 1) ? 2 : 1;
    } else {
        be.jacobian_type = jacobian;
    }

    // Start off by assuming a successful burn.

    state.success = true;

    // Initialize the integration time.

    be.t = 0.0_rt;
    be.tout = dt;

    // Initialize ydot to zero for Strang burn.

    for (int n = 0; n < SVAR; ++n) {
        state.ydot_a[n] = 0;
    }

    // We assume that (rho, T) coming in are valid, do an EOS call
    // to fill the rest of the thermodynamic variables.

    eos(eos_input_rt, state);

    // set the scaling for energy if we integrate it dimensionlessly
    state.e_scale = state.e;

    if (scale_system) {
        // the absolute tol for energy needs to reflect the scaled
        // energy the integrator sees
        be.atol_enuc /= state.e_scale;
    }

    // Fill in the initial integration state.

    burn_to_integrator(state, be);

    // Save the initial composition, temperature, and energy for our later diagnostics.

#ifndef AMREX_USE_GPU
    amrex::Real xn_in[NumSpec];
    for (int n = 0; n < NumSpec; ++n) {
        xn_in[n] = state.xn[n];
    }
    const amrex::Real T_in = state.T;
#endif
    const amrex::Real e_in = state.e;

    // Call the integration routine.

    int istate = be_integrator(state, be);
    state.error_code = istate;

    // Copy the integration data back to the burn state.

    integrator_to_burn(be, state);

#ifdef NSE
    // compute the temperature based on the energy release -- we need
    // this in case we failed in our burn here because we entered NSE

#ifdef AUX_THERMO
    // need to sync the auxiliary data up with the new mass fractions
    set_aux_comp_from_X(state);
#endif
    if (call_eos_in_rhs) {
        eos(eos_input_re, state);
    }
#endif

    // Subtract off the initial energy if the application codes expect
    // to get back only the generated energy during the burn.
    if (integrator_rp::subtract_internal_energy) {
        state.e -= e_in;
    }

    // Normalize the final abundances.

    if (! integrator_rp::use_number_densities) {
        normalize_abundances_burn(state);
    }

    // Get the number of RHS and Jacobian evaluations.

    state.n_rhs = be.n_rhs;
    state.n_jac = be.n_jac;
    state.n_step = be.n_step;

    // BE does not always fail even though it can lead to unphysical states.
    // Add some checks that indicate a burn fail even if VODE thinks the
    // integration was successful.

    if (istate != IERR_SUCCESS) {
        state.success = false;
    }

    for (int n = 1; n <= NumSpec; ++n) {
        if (be.y(n) < -species_failure_tolerance) {
            state.success = false;
        }

        // Don't enforce a max if we are evolving number densities

        if (! integrator_rp::use_number_densities) {
            if (be.y(n) > 1.0_rt + species_failure_tolerance) {
                state.success = false;
            }
        }
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
        if (istate != IERR_ENTERED_NSE) {
#ifndef AMREX_USE_GPU
            std::cout << amrex::Font::Bold << amrex::FGColor::Red << "[ERROR] integration failed in net" << amrex::ResetDisplay << std::endl;
            std::cout << "istate = " << istate << std::endl;
            if (istate == IERR_SUCCESS) {
                std::cout << "  BE exited successfully, but a check on the data values failed" << std::endl;
            }
            std::cout << "zone = (" << state.i << ", " << state.j << ", " << state.k << ")" << std::endl;
            std::cout << "time = " << be.t << std::endl;
            std::cout << "dt = " << std::setprecision(16) << dt << std::endl;
            std::cout << "temp start = " << std::setprecision(16) << T_in << std::endl;
            std::cout << "xn start = ";
            for (const auto X: xn_in) {
                std::cout << std::setprecision(16) << X << " ";
            }
            std::cout << std::endl;
            std::cout << "dens current = " << std::setprecision(16) << state.rho << std::endl;
            std::cout << "temp current = " << std::setprecision(16) << state.T << std::endl;
            std::cout << "xn current = ";
            for (const auto X: state.xn) {
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
