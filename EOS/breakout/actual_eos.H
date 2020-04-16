#ifndef _actual_eos_H_
#define _actual_eos_H_

#include <extern_parameters.H>
#include <fundamental_constants.H>
#include <network.H>
#include <actual_eos_data.H>

const std::string eos_name = "breakout";

inline
void actual_eos_init ()
{

    // constant ratio of specific heats
    if (eos_gamma > 0.e0_rt) {
        gamma_const = eos_gamma;
    } else {
        gamma_const = 5.0_rt / 3.0_rt;
    }

}

AMREX_GPU_HOST_DEVICE inline
bool is_input_valid(eos_input_t input) {

  bool valid = true;

  if (input == eos_input_rh ||
      input == eos_input_tp ||
      input == eos_input_ps ||
      input == eos_input_ph ||
      input == eos_input_th) {
    valid = false;
  }

  return valid;
}

AMREX_GPU_HOST_DEVICE inline
void actual_eos (eos_input_t input, eos_t& state)
{

    const Real R = C::k_B * C::n_A;

    Real poverrho;

    // Calculate mu. This is the only difference between
    // this EOS and gamma_law.
    state.mu = 1.0_rt / state.aux[1];

    switch (input) {

    case eos_input_rt:

        // dens, temp and xmass are inputs
        state.cv = R / (state.mu * (gamma_const - 1.0_rt)) ;
        state.e = state.cv * state.T;
        state.p = (gamma_const - 1.0_rt) * state.rho * state.e;
        state.gam1 = gamma_const;

        break;

    case eos_input_rh:

        // dens, enthalpy, and xmass are inputs
#ifndef AMREX_USE_GPU
        amrex::Error("EOS: eos_input_rh is not supported in this EOS.");
#endif

        break;

    case eos_input_tp:

        // temp, pres, and xmass are inputs
#ifndef AMREX_USE_GPU
        amrex::Error("EOS: eos_input_tp is not supported in this EOS.");
#endif

        break;

    case eos_input_rp:

        // dens, pres, and xmass are inputs

        poverrho = state.p / state.rho;
        state.T = poverrho * state.mu * (1.0_rt / R);
        state.e = poverrho * (1.0_rt / (gamma_const - 1.0_rt));
        state.gam1 = gamma_const;

        break;

    case eos_input_re:

        // dens, energy, and xmass are inputs

        poverrho = (gamma_const - 1.0_rt) * state.e;

        state.p = poverrho * state.rho;
        state.T = poverrho * state.mu * (1.0_rt / R);
        state.gam1 = gamma_const;

        // sound speed
        state.cs = sqrt(gamma_const * poverrho);

        state.dpdr_e = poverrho;
        state.dpde = (gamma_const - 1.0_rt) * state.rho;

        // Try to avoid the expensive log function.  Since we don't need entropy 
        // in hydro solver, set it to an invalid but "nice" value for the plotfile.
        state.s = 1.0_rt;

        break;

    case eos_input_ps:

        // pressure entropy, and xmass are inputs

#ifndef AMREX_USE_GPU
        amrex::Error("EOS: eos_input_ps is not supported in this EOS.");
#endif

        break;

    case eos_input_ph:

        // pressure, enthalpy and xmass are inputs
#ifndef AMREX_USE_GPU
        amrex::Error("EOS: eos_input_ph is not supported in this EOS.");
#endif

        break;

    case eos_input_th:

        // temperature, enthalpy and xmass are inputs

        // This system is underconstrained.
#ifndef AMREX_USE_GPU
        amrex::Error("EOS: eos_input_th is not a valid input for the gamma law EOS.");
#endif

        break;

    default:

#ifndef AMREX_USE_GPU
        amrex::Error("EOS: invalid input.");
#endif

        break;

    }

}



inline
void actual_eos_finalize ()
{
}

#endif
