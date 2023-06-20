#ifndef actual_integrator_H
#define actual_integrator_H

#include <vode_type.H>
#include <vode_util.H>
#include <core_integrator.H>

template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void actual_integrator (BurnT& state, Real dt)
{

    constexpr int int_neqs = integrator_neqs<BurnT>();

    dvode_t<int_neqs> int_state{};

    core_integrator(state, int_state, dt);
}

#endif
