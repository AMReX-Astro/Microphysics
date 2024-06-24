#ifndef actual_integrator_H
#define actual_integrator_H

#include <network.H>
#include <burn_type.H>

#include <integrator_setup_sdc.H>

#include <be_type.H>
#include <be_integrator.H>

template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_integrator (BurnT& state, const amrex::Real dt, bool is_retry=false)
{
    constexpr int int_neqs = integrator_neqs<BurnT>();

    auto be_state = integrator_setup<BurnT, be_t<int_neqs>>(state, dt, is_retry);
    auto state_save = integrator_backup(state);

    // Call the integration routine.

    int istate = be_integrator(state, be_state);
    state.error_code = istate;

    integrator_cleanup(be_state, state, istate, state_save, dt);

}

#endif
