#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#ifdef SIMPLIFIED_SDC
#include <actual_integrator_simplified_sdc.H>
#else
#include <actual_integrator.H>
#endif

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void integrator (burn_t& state, Real dt)
{
    actual_integrator(state, dt);
}

#endif
