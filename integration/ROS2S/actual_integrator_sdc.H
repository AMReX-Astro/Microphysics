#ifndef actual_integrator_H
#define actual_integrator_H

#include <AMReX_REAL.H>

#include <network.H>
#include <burn_type.H>

#include <integrator_data.H>
#include <integrator_setup_sdc.H>
#include <ros2s_integrator.H>
#include <ros2s_type.H>

template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_integrator (BurnT& state, amrex::Real dt, bool is_retry=false)
{
    constexpr int int_neqs = integrator_neqs<BurnT>();

    auto ros2s_state = integrator_setup<BurnT, ros2s_t<int_neqs>>(state, dt, is_retry);
    auto state_save = integrator_backup(state);

    int istate = ros2s_integrator(state, ros2s_state);

    integrator_cleanup(ros2s_state, state, istate, state_save, dt);

    state.rhs_walltime = ros2s_state.rhs_walltime;
    state.jacobian_walltime = ros2s_state.jacobian_walltime;
    state.lu_walltime = ros2s_state.lu_walltime;
    state.linear_solve_walltime = ros2s_state.linear_solve_walltime;
}

#endif
