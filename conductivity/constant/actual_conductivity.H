#ifndef _actual_conductivity_H_
#define _actual_conductivity_H_

#include <cmath>
#include "eos_type.H"
#include "network.H"
#include "actual_conductivity_data.H"


AMREX_FORCE_INLINE
void
actual_conductivity_init()
{
  const_cond = const_conductivity;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
actual_conductivity(eos_t state) {
  eos_state.conductivity = const_cond;
}
