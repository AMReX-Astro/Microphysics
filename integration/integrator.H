#ifndef _integrator_H_
#define _integrator_H_

#include <actual_integrator.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void integrator (burn_t& state, Real dt)
{
    actual_integrator(state, dt);
}

#endif
