#ifndef SDC_INTEGRATOR_H
#define SDC_INTEGRATOR_H

#include <vode_integrator.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void integrator (burn_t& state, Real dt);
{
  vode_integrator(state, dt);
}

#endif
