#ifndef EOS_COMPOSITION_H
#define EOS_COMPOSITION_H

#include <AMReX_BLFort.H>
#include <network.H>
#include <eos_type.H>

#ifdef AUX_THERMO
#include <actual_network.H>
#endif

using namespace amrex;

#ifdef AUX_THERMO
using namespace AuxZero;
#endif

struct eos_xderivs_t {
  Real dedX[NumSpec];
  Real dpdX[NumSpec];
  Real dhdX[NumSpec];
};

#ifdef AUX_THERMO
///
/// The auxiliary state provides an alternate description to the composition,
/// in terms of Ye, abar, and binding energy / nucleon
///
template <class state_t> 
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void set_aux_comp_from_X(state_t& state) {

  using namespace AuxZero;

  state.aux[iye] = 0.0_rt;
  state.aux[iabar] = 0.0_rt;
  state.aux[ibea] = 0.0_rt;

  for (int n = 0; n < NumSpec; n++) {
    state.aux[iye] += state.xn[n] * zion[n] * aion_inv[n];
    state.aux[iabar] += state.xn[n] * aion_inv[n];
    state.aux[ibea] += state.xn[n] * network::bion(n+1) * aion_inv[n];
  }
  state.aux[iabar] = 1.0_rt/state.aux[iabar];

}
#endif

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void composition (T& state)
{
  // Calculate abar, the mean nucleon number,
  // zbar, the mean proton number,
  // mu, the mean molecular weight,
  // mu_e, the mean number of nucleons per electron, and
  // y_e, the electron fraction.

#ifdef AUX_THERMO

  state.mu_e = 1.0 / state.aux[iye];
  state.y_e = state.aux[iye];
  state.abar = state.aux[iabar];
  state.zbar = state.abar * state.y_e;

#else

  Real sum = 0;
  for (int n = 0; n < NumSpec; n++) {
    sum += state.xn[n] * zion[n] * aion_inv[n];
  }

  state.mu_e = 1.0 / sum;
  state.y_e = sum;

  sum = 0;

  for (int n = 0; n < NumSpec; n++) {
    sum += state.xn[n] * aion_inv[n];
  }

  state.abar = 1.0 / sum;
  state.zbar = state.abar * state.y_e;
#endif
}

#ifdef EXTRA_THERMO
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
eos_xderivs_t composition_derivatives (const T& state)
{
  eos_xderivs_t state_xderivs;

  for (int n = 0; n < NumSpec; n++) {
      if constexpr (has_pressure<T>::value) {
          state_xderivs.dpdX[n] =
              state.dpdA * (state.abar * aion_inv[n]) * (aion[n] - state.abar) +
              state.dpdZ * (state.abar * aion_inv[n]) * (zion[n] - state.zbar);
      }

      if constexpr (has_energy<T>::value) {
          state_xderivs.dedX[n] =
              state.dedA * (state.abar * aion_inv[n]) * (aion[n] - state.abar) +
              state.dedZ * (state.abar * aion_inv[n]) * (zion[n] - state.zbar);
      }

      if constexpr (has_pressure<T>::value) {
          if (state.dpdr != 0.0) {

              state_xderivs.dhdX[n] = state_xderivs.dedX[n]
                  + (state.p / (state.rho*state.rho) - state.dedr)
                  * state_xderivs.dpdX[n] / state.dpdr;

          }
      }
  }

  return state_xderivs;
}
#endif


#endif
