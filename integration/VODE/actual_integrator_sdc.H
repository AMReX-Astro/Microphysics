#ifndef actual_integrator_H
#define actual_integrator_H

// Common variables and routines for burners
// that use VODE for their integration.

#include <network.H>
#include <burn_type.H>

#include <integrator_setup_sdc.H>

#include <vode_type.H>
#include <vode_dvode.H>

using namespace integrator_rp;

template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_integrator (BurnT& state, amrex::Real dt, bool is_retry=false)
{

    constexpr int int_neqs = integrator_neqs<BurnT>();

    auto vode_state = integrator_setup<BurnT, dvode_t<int_neqs>>(state, dt, is_retry);
    auto state_save = integrator_backup(state);

    // Call the integration routine.

    auto istate = dvode(state, vode_state);
    state.error_code = istate;

    integrator_cleanup(vode_state, state, istate, state_save, dt);


}

#endif
