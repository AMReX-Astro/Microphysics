#ifndef _actual_conductivity_H_
#define _actual_conductivity_H_

#include <cmath>
#include "eos_type.H"
#include "network.H"

const std::string cond_name = "powerlaw";

AMREX_FORCE_INLINE
void
actual_conductivity_init()
{
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
actual_conductivity(eos_t& state) {
  state.conductivity = cond_coeff * std::pow(state.T, cond_exponent);
}
#endif
