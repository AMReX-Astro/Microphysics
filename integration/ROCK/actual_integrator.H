#ifndef actual_integrator_H
#define actual_integrator_H

#include <network.H>
#include <burn_type.H>

#include <integrator_data.H>
#include <integrator_setup_strang.H>

#include <rock_type.H>
#include <rock.H>

template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_integrator (BurnT& state, amrex::Real dt, bool is_retry=false)
{
    constexpr int int_neqs = integrator_neqs<BurnT>();

    auto rock_state = integrator_setup<BurnT, rock_t<int_neqs>>(state, dt, is_retry);

    auto state_save = integrator_backup(state);

    auto istate = rock(state, rock_state);

    integrator_cleanup(rock_state, state, istate, state_save, dt);

}

#endif
