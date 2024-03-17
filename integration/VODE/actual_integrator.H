#ifndef actual_integrator_H
#define actual_integrator_H

// Common variables and routines for burners
// that use VODE for their integration.

#include <iomanip>

#include <network.H>
#include <burn_type.H>
#include <eos_type.H>
#include <eos.H>
#include <extern_parameters.H>
#include <integrator_data.H>

#include <integrator_setup.H>

#include <vode_type.H>
#include <vode_dvode.H>

template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_integrator (BurnT& state, amrex::Real dt, bool is_retry=false)
{

    constexpr int int_neqs = integrator_neqs<BurnT>();

    integrator_setup<BurnT, dvode_t<int_neqs>>(state, dt, is_retry,
                                               dvode<BurnT, dvode_t<int_neqs>>);

}

#endif
