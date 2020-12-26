#ifndef _eos_override_H_
#define _eos_override_H_

#include <eos_type.H>

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void eos_override (T& /*state*/) {}

#endif

