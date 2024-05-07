#ifndef actual_integrator_H
#define actual_integrator_H

#include <network.H>
#include <burn_type.H>

#include <integrator_data.H>
#include <integrator_setup.H>

#include <be_type.H>
#include <be_integrator.H>

template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_integrator (BurnT& state, const amrex::Real dt, bool is_retry=false)
{

    constexpr int int_neqs = integrator_neqs<BurnT>();

    integrator_setup<BurnT, be_t<int_neqs>>(state, dt, is_retry,
                                            be_integrator<BurnT, be_t<int_neqs>>);

}

#endif
