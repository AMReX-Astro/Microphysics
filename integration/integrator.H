#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#ifdef SIMPLIFIED_SDC
#include <actual_integrator_simplified_sdc.H>
#else
#include <actual_integrator.H>
#endif

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void integrator (burn_t& state, Real dt)
{
    actual_integrator(state, dt);

#ifndef NSE_TABLE
    // optionally abort here if we don't have retries in the driver.
    // Note, for NSE we get a second chance if the state entered NSE,
    // in burner()
    if (!state.success && abort_on_failure) {
#ifndef AMREX_USE_GPU
        std::cout << "failure: \n";
        std::cout << state << std::endl;
        amrex::Error("integration failed in burner");
#endif
    }
#endif

}

#endif
