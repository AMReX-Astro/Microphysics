#ifndef actual_integrator_H
#define actual_integrator_H

#include <AMReX_Print.H>

#include <iomanip>

#include <network.H>
#include <burn_type.H>
#include <eos_type.H>
#include <eos.H>
#include <extern_parameters.H>
#include <rkc_type.H>
#include <rkc.H>

template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_integrator (BurnT& state, amrex::Real dt, bool is_retry=false)
{
    constexpr int int_neqs = integrator_neqs<BurnT>();

    rkc_t<int_neqs> rkc_state{};

    // Set the tolerances.

    if (!is_retry) {
        rkc_state.atol_spec = atol_spec;  // mass fractions
        rkc_state.atol_enuc = atol_enuc;  // energy generated

        rkc_state.rtol_spec = rtol_spec;  // mass fractions
        rkc_state.rtol_enuc = rtol_enuc;  // energy generated
    } else {
        rkc_state.atol_spec = retry_atol_spec; // mass fractions
        rkc_state.atol_enuc = retry_atol_enuc; // energy generated

        rkc_state.rtol_spec = retry_rtol_spec; // mass fractions
        rkc_state.rtol_enuc = retry_rtol_enuc; // energy generated
    }


    // Start off by assuming a successful burn.

    state.success = true;

    // Initialize the integration time.

    rkc_state.t = 0.0_rt;
    rkc_state.tout = dt;

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
        rkc_state.atol_enuc /= state.e_scale;
    }

    // Fill in the initial integration state.

    burn_to_integrator(state, rkc_state);

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

    int istate = rkc(state, rkc_state);
    state.error_code = istate;

    // Copy the integration data back to the burn state.

    integrator_to_burn(rkc_state, state);

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

    // Subtract off the initial energy (the application codes expect
    // to get back only the generated energy during the burn).
    // Don't subtract it for primordial chem
    if (subtract_internal_energy) {
        state.e -= e_in;
    }

    // Normalize the final abundances.
    // Don't normalize for primordial chem
    if (!use_number_densities) {
        normalize_abundances_burn(state);
    }

    // Get the number of RHS and Jacobian evaluations.

    state.n_rhs = rkc_state.nfe;
    state.n_jac = 0;
    state.n_step = rkc_state.nsteps;

    if (istate != IERR_SUCCESS) {
        state.success = false;
    }

    // We may not fail even though we generate unphysical states.
    // Add some checks that indicate a burn fail even if RKC thinks the
    // integration was successful.

    for (int n = 1; n <= NumSpec; ++n) {
        if (rkc_state.y(n) < -species_failure_tolerance) {
            state.success = false;
        }

        // Don't enforce the condition below
        // for primordial chem
        if (!use_number_densities) {
            if (rkc_state.y(n) > 1.0_rt + species_failure_tolerance) {
                state.success = false;
            }
        }
    }

#ifndef AMREX_USE_GPU
    if (burner_verbose) {
        // Print out some integration statistics, if desired.
        std::cout <<  "integration summary: " << std::endl;
        std::cout <<  "dens: " << state.rho << " temp: " << state.T << std::endl;
        std::cout << " energy released: " << state.e << std::endl;
        std::cout <<  "number of steps taken: " << rkc_state.nsteps << std::endl;
        std::cout <<  "number of f evaluations: " << rkc_state.nfe << std::endl;
    }
#endif

    // If we failed, print out the current state of the integration.

    if (!state.success) {
        if (istate != IERR_ENTERED_NSE) {
#ifndef AMREX_USE_GPU
            std::cout << amrex::Font::Bold << amrex::FGColor::Red << "[ERROR] integration failed in net" << amrex::ResetDisplay << std::endl;
            std::cout << "istate = " << istate << std::endl;
            if (istate == IERR_SUCCESS) {
                std::cout << "  RKC exited successfully, but a check on the data values failed" << std::endl;
            }
            std::cout << "zone = (" << state.i << ", " << state.j << ", " << state.k << ")" << std::endl;
            std::cout << "time = " << rkc_state.t << std::endl;
            std::cout << "dt = " << std::setprecision(16) << dt << std::endl;
            std::cout << "temp start = " << std::setprecision(16) << T_in << std::endl;
            std::cout << "xn start = ";
            for (auto X : xn_in) {
                std::cout << std::setprecision(16) << X << " ";
            }
            std::cout << std::endl;
            std::cout << "dens current = " << std::setprecision(16) << state.rho << std::endl;
            std::cout << "temp current = " << std::setprecision(16) << state.T << std::endl;
            std::cout << "xn current = ";
            for (auto X : state.xn) {
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
