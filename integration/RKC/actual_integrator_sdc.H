#ifndef actual_integrator_H
#define actual_integrator_H

#include <network.H>
#include <burn_type.H>

#include <integrator_setup_sdc.H>

#include <rkc_type.H>
#include <rkc.H>

using namespace integrator_rp;

template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_integrator (BurnT& state, amrex::Real dt, bool is_retry=false)
{

    constexpr int int_neqs = integrator_neqs<BurnT>();

    auto rkc_state = integrator_setup<BurnT, rkc_t<int_neqs>>(state, dt, is_retry);
    auto state_save = integrator_backup(state);

    // Call the integration routine.

    int istate = rkc(state, rkc_state);
    state.error_code = istate;

    integrator_cleanup(rkc_state, state, istate, state_save, dt);

}
#endif
