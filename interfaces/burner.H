#ifndef _burner_H_
#define _burner_H_

#include <burn_type.H>
#ifdef SIMPLIFIED_SDC
#include <integrator_sdc.H>
#else
#include <integrator.H>
#endif

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void burner (burn_t& state, Real dt)
{
    integrator(state, dt);
}

#endif
