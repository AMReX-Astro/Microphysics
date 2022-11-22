#ifndef actual_integrator_H
#define actual_integrator_H

#include <be_type.H>
#include <be_integrator.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void actual_integrator (burn_t& state, const Real dt)
{

    // start off by assuming a successful burn

    state.success = true;

    be_t be;

    be.t = 0.0;
    be.tout = dt;

    be.atol_spec = atol_spec;
    be.rtol_spec = rtol_spec;

    be.atol_enuc = atol_enuc;
    be.rtol_enuc = rtol_enuc;

    // We assume that (rho, T) coming in are valid, do an EOS call
    // to fill the rest of the thermodynamic variables.

    eos(eos_input_rt, state);

    // Fill in the initial integration state.

    burn_to_integrator(state, be);

#ifndef AMREX_USE_GPU
    Real T_in = state.T;
    Real xn_in[NumSpec];
    for (int n = 0; n < NumSpec; ++n) {
        xn_in[n] = state.xn[n];
    }
#endif
    Real e_in = state.e;

    // Call the integration routine.

    int ierr = be_integrator(state, be);

    // Copy the integration data back to the burn state.

    integrator_to_burn(be, state);

    // Subtract off the initial energy if the application codes expect
    // to get back only the generated energy during the burn.
    if (subtract_internal_energy) {
        state.e -= e_in;
    }

    // Normalize the final abundances.

    normalize_abundances_burn(state);

    state.n_rhs = be.n_rhs;
    state.n_step = be.n_step;

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
