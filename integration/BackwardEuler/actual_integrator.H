#ifndef actual_integrator_H
#define actual_integrator_H

#include <be_integrator.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void actual_integrator (burn_t& state, const Real dt)
{
    be_integrator(state, dt);
}

#endif
