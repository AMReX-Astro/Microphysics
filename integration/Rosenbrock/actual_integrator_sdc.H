#ifndef actual_integrator_H
#define actual_integrator_H

#include <AMReX_REAL.H>

#include <network.H>
#include <burn_type.H>

#include <integrator_data.H>
#include <integrator_setup_sdc.H>
#include <rosenbrock_integrator.H>
#include <rosenbrock_type.H>

template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_integrator (BurnT& state, amrex::Real dt, bool is_retry=false)
{
    constexpr int int_neqs = integrator_neqs<BurnT>();

    if (integrator_rp::rosenbrock_tableau == 0) {
        auto rodas5p_state = integrator_setup<BurnT, rodas5p_t<int_neqs>>(state, dt, is_retry);
        auto state_save = integrator_backup(state);

        int istate = rodas5p_integrator(state, rodas5p_state);

        integrator_cleanup(rodas5p_state, state, istate, state_save, dt);
        return;
    }

    auto rosenbrock_state = integrator_setup<BurnT, rosenbrock_t<int_neqs>>(state, dt, is_retry);
    auto state_save = integrator_backup(state);

    int istate = rosenbrock_integrator(state, rosenbrock_state);

    integrator_cleanup(rosenbrock_state, state, istate, state_save, dt);
}

#endif
