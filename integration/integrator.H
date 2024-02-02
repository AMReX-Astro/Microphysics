#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#ifdef SDC
#include <actual_integrator_simplified_sdc.H>
#else
#include <actual_integrator.H>
#endif

template <typename BurnT, bool enable_retry>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void integrator_wrapper (BurnT& state, Real dt)
{

    if constexpr (enable_retry) {
        burn_t old_state{state};

        actual_integrator(state, dt);

        if (!state.success) {
            state = old_state;
            const bool is_retry = true;
            actual_integrator(state, dt, is_retry);
        }
    } else {
        actual_integrator(state, dt);
    }

}


template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void integrator (BurnT& state, Real dt)
{

    if (integrator_rp::use_burn_retry) {
        constexpr bool enable_retry{true};
        integrator_wrapper<BurnT, enable_retry>(state, dt);
    } else {
        constexpr bool enable_retry{false};
        integrator_wrapper<BurnT, enable_retry>(state, dt);

    }
}

#endif
