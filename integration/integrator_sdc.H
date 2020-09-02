#ifndef SDC_INTEGRATOR_H
#define SDC_INTEGRATOR_H

#include <vode_integrator.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void integrator (burn_t& state_in, burn_t& state_out, Real dt, Real time)
{
  vode_integrator(state_in, state_out, dt, time)
}

#endif
