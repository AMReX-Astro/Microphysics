#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#ifdef SDC
#include <actual_integrator_simplified_sdc.H>
#else
#include <actual_integrator.H>
#endif

template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void integrator (BurnT& state, amrex::Real dt)
{

    // create a backup in case we do a retry

    burn_t old_state{state};

    actual_integrator(state, dt);

    if (integrator_rp::use_burn_retry) {
        if (!state.success) {
            state = old_state;
            const bool is_retry = true;
            actual_integrator(state, dt, is_retry);
        }
    }
}

#endif
