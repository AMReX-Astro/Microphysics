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



AMREX_GPU_HOST_DEVICE inline
void actual_eos (eos_input_t input, eos_t& state)
{

    switch (input) {

    case eos_input_rt:

        state.cv = eos_const_c_v * std::pow(state.rho, eos_c_v_exp_m) * std::pow(state.T, -eos_c_v_exp_n);
        state.e  = eos_const_c_v * std::pow(state.rho, eos_c_v_exp_m) * std::pow(state.T, 1 - eos_c_v_exp_n) / (1 - eos_c_v_exp_n);

        break;

    case eos_input_re:

        state.T  = std::pow((1 - eos_c_v_exp_n) * state.e * std::pow(state.rho, -eos_c_v_exp_m) / eos_const_c_v, 1 / (1 - eos_c_v_exp_n));
        state.cv = eos_const_c_v * std::pow(state.rho, eos_c_v_exp_m) * std::pow(state.T, -eos_c_v_exp_n);

        break;

    default:

#ifndef AMREX_USE_GPU
        amrex::Error("EOS: invalid input.");
#endif

        break;

    }

    // Set some data to nonsense values so that things intentionally go wrong
    // if this EOS is somehow used for hydro.

    state.p    = -1.e0_rt;
    state.gam1 = -1.e0_rt;
    state.cs   = -1.e0_rt;
    state.s    = -1.e0_rt;
    state.h    = -1.e0_rt;

    // Give dpdr a value for the purposes of the composition_derivatives routine.

    state.dpdr = 0.e0_rt;

}

#endif
