#ifndef ACTUAL_CONDUCTIVITY_H
#define ACTUAL_CONDUCTIVITY_H

#include <string_view>
#include <cmath>
#include <eos_type.H>
#include <network.H>
#include <extern_parameters.H>

using namespace conductivity_rp;

constexpr std::string_view cond_name = "constant_opacity";

AMREX_FORCE_INLINE
void
actual_conductivity_init()
{
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
actual_conductivity (T& state)
{
  state.conductivity = (16*sigma_SB*state.T*state.T*state.T)/(3*const_opacity*state.rho);
}

#endif
