#ifndef SDC_INTEGRATOR_H
#define SDC_INTEGRATOR_H

#include <actual_integrator_simplified_sdc.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void integrator (burn_t& state, Real dt)
{
    actual_integrator(state, dt);
}

#endif
