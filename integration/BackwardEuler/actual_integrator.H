#ifndef actual_integrator_H
#define actual_integrator_H

#include <be_type.H>
#include <be_integrator.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void actual_integrator (burn_t& state, const Real dt)
{

    be_t be;

    // Set the tolerances.

    be.atol_spec = atol_spec;  // mass fractions
    be.atol_enuc = atol_enuc;  // energy generated

    be.rtol_spec = rtol_spec;  // mass fractions
    be.rtol_enuc = rtol_enuc;  // energy generated

    // set the Jacobian type
    be.jacobian_type = jacobian;

    // Start off by assuming a successful burn.

    state.success = true;

    // Initialize the integration time.

    be.t = 0.0_rt;
    be.tout = dt;


    // We assume that (rho, T) coming in are valid, do an EOS call
    // to fill the rest of the thermodynamic variables.

    eos(eos_input_rt, state);

    // Fill in the initial integration state.

    burn_to_integrator(state, be);

    // Save the initial composition, temperature, and energy for our later diagnostics.

#ifndef AMREX_USE_GPU
    Real xn_in[NumSpec];
    for (int n = 0; n < NumSpec; ++n) {
        xn_in[n] = state.xn[n];
    }
    Real T_in = state.T;
#endif
    Real e_in = state.e;

    // Call the integration routine.

    int istate = be_integrator(state, be);

    // Copy the integration data back to the burn state.

    integrator_to_burn(be, state);

    // Subtract off the initial energy if the application codes expect
    // to get back only the generated energy during the burn.
    if (subtract_internal_energy) {
        state.e -= e_in;
    }

    // Normalize the final abundances.

    if (!use_number_densities) {
        normalize_abundances_burn(state);
    }

    // Get the number of RHS and Jacobian evaluations.

    state.n_rhs = be.n_rhs;
    state.n_jac = be.n_jac;
    state.n_step = be.n_step;

    if (istate < 0) {
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
        std::cout << Font::Bold << FGColor::Red << "[ERROR] integration failed in net" << ResetDisplay << std::endl;
        std::cout << "istate = " << istate << std::endl;
        std::cout << "zone = (" << state.i << ", " << state.j << ", " << state.k << ")" << std::endl;
        std::cout << "time = " << be.t << std::endl;
        std::cout << "dt = " << std::setprecision(16) << dt << std::endl;
        std::cout << "temp start = " << std::setprecision(16) << T_in << std::endl;
        std::cout << "xn start = ";
        for (int n = 0; n < NumSpec; ++n) {
            std::cout << std::setprecision(16) << xn_in[n] << " ";
        }
        std::cout << std::endl;
        std::cout << "dens current = " << std::setprecision(16) << state.rho << std::endl;
        std::cout << "temp current = " << std::setprecision(16) << state.T << std::endl;
        std::cout << "xn current = ";
        for (int n = 0; n < NumSpec; ++n) {
            std::cout << std::setprecision(16) << state.xn[n] << " ";
        }
        std::cout << std::endl;
        std::cout << "energy generated = " << state.e << std::endl;
#endif
    }

}

#endif
