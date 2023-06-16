#ifndef CONDUCTIVITY_H
#define CONDUCTIVITY_H

#include <eos_type.H>
#include <actual_conductivity.H>

using namespace amrex;

AMREX_FORCE_INLINE
void conductivity_init() {
  actual_conductivity_init();
}


template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void conductivity (T& state)
{
  actual_conductivity(state);
}
#endif
