#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#ifdef SDC
#include <actual_integrator_simplified_sdc.H>
#else
#include <actual_integrator.H>
#endif

template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void integrator (BurnT& state, Real dt)
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
}

#endif
