#ifndef _actual_eos_H_
#define _actual_eos_H_

// This is an artificial equation of state used primarily for radiation tests.
//
// It is defined by the relationship:
// c_v = K * rho**m * T**(-n)
// where K, m, and n are user-defined constant parameters.
//
// Ignoring the integration constant, we thus have the relationships:
// e = K * rho**m * T**(1-n) / (1 - n)
// T = ((1 - n) * e * rho**(-m) / K)**(1/(1-n))
//
// Consequently the only input modes supported are eos_input_rt and eos_input_re.
// Pressure and Gamma_1 are not defined, so this EOS cannot be used for hydro.

#include <extern_parameters.H>
#include <eos_type.H>
#include <cmath>

const std::string eos_name = "rad_power_law";

inline
void actual_eos_init ()
{

    if (eos_const_c_v <= 0.e0_rt) {
        amrex::Error("eos_const_c_v must be > 0");
    }

    if (eos_c_v_exp_n == 1.0e0_rt) {
        amrex::Error("eos_c_v_exp_n == 1 is unsupported");
    }

}



inline
void actual_eos_finalize ()
{
}



template <typename I>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
bool is_input_valid (I input)
{
  static_assert(std::is_same<I, eos_input_t>::value, "input must be an eos_input_t");

  bool valid = true;

  if (input == eos_input_rh ||
      input == eos_input_tp ||
      input == eos_input_rp ||
      input == eos_input_ps ||
      input == eos_input_ph ||
      input == eos_input_th) {
    valid = false;
  }

  return valid;
}


template <typename I, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_eos (I input, T& state)
{
    static_assert(std::is_same<I, eos_input_t>::value, "input must be an eos_input_t");

    switch (input) {

    case eos_input_rt:

        if constexpr (has_energy<T>::value) {
            state.cv = eos_const_c_v * std::pow(state.rho, eos_c_v_exp_m) * std::pow(state.T, -eos_c_v_exp_n);
            state.e  = eos_const_c_v * std::pow(state.rho, eos_c_v_exp_m) * std::pow(state.T, 1 - eos_c_v_exp_n) / (1 - eos_c_v_exp_n);
        }

        break;

    case eos_input_re:

        if constexpr (has_energy<T>::value) {
            state.T  = std::pow((1 - eos_c_v_exp_n) * state.e * std::pow(state.rho, -eos_c_v_exp_m) / eos_const_c_v, 1.0_rt / (1.0_rt - eos_c_v_exp_n));
            state.cv = eos_const_c_v * std::pow(state.rho, eos_c_v_exp_m) * std::pow(state.T, -eos_c_v_exp_n);
        }

        break;

    default:

#ifndef AMREX_USE_GPU
        amrex::Error("EOS: invalid input.");
#endif

        break;

    }

    // Set some data to nonsense values so that things intentionally go wrong
    // if this EOS is somehow used for hydro.

    if constexpr (has_pressure<T>::value) {
        state.p    = -1.e0_rt;
    }
    if constexpr (has_pressure<T>::value && has_energy<T>::value) {
        state.gam1 = -1.e0_rt;
        state.cs   = -1.e0_rt;
    }
    if constexpr (has_entropy<T>::value) {
        state.s    = -1.e0_rt;
    }
    if constexpr (has_enthalpy<T>::value) {
        state.h    = -1.e0_rt;
    }

    // Give dpdr a value for the purposes of the composition_derivatives routine.

    if constexpr (has_pressure<T>::value) {
        state.dpdr = 0.e0_rt;
    }
}

#endif
