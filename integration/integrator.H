#ifndef _integrator_H_
#define _integrator_H_

#ifdef SIMPLIFIED_SDC
#include <actual_integrator_simplified_sdc.H>
#else
#include <actual_integrator.H>
#endif

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void integrator (burn_t& state, Real dt)
{

    burn_t old_state;
    if (use_burn_retry) {
        old_state = state;
    }

    actual_integrator(state, dt);

    if (!state.success && use_burn_retry) {
        state = old_state;
        const bool is_retry = true;
        actual_integrator(state, dt, is_retry);
    }

#ifndef NSE_TABLE
    if (!state.success && abort_on_failure) {
        // optionally abort here if we don't have retries in the driver.
        // Note, for NSE we get a second chance if the state entered NSE,
        // in burner()

#ifndef AMREX_USE_GPU
        std::cout << "failure: \n";
        std::cout << state << std::endl;
        amrex::Error("integration failed in burner");
#endif
    }
#endif

}

#endif
