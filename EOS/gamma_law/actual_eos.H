#ifndef _actual_eos_H_
#define _actual_eos_H_

// This is a constant gamma equation of state, using an ideal gas.
//
// This a simplified version of the more general eos_gamma_general.

#include <actual_eos_data.H>
#include <extern_parameters.H>
#include <fundamental_constants.H>
#include <network.H>
#include <eos_type.H>
#include <cmath>

const std::string eos_name = "gamma_law";

inline
void actual_eos_init ()
{

    // constant ratio of specific heats
    if (eos_gamma > 0.e0_rt)
    {
        gamma_const = eos_gamma;
    }
    else
    {
        amrex::Error("gamma_const cannot be < 0");
    }

    assume_neutral = eos_assume_neutral;

}

inline
void actual_eos_finalize ()
{
}

AMREX_GPU_HOST_DEVICE inline
void actual_eos (eos_input_t input, eos_t& state)
{

    const Real R = k_B * n_A;

    // Calculate mu.

    if (assume_neutral)
    {
        state.mu = state.abar;
    }
    else
    {
        Real sum = 0.0_rt;
        for (int n = 0; n < NumSpec; ++n) {
            sum = sum + (1.0_rt + zion[n]) * state.xn[n] * aion_inv[n];
        }
        state.mu = 1.0_rt / sum;
    }

    Real poverrho;

    switch (input) {

    case eos_input_rt:

        // dens, temp and xmass are inputs
        state.cv = R / (state.mu * (gamma_const - 1.0_rt));
        state.e = state.cv * state.T;
        state.p = (gamma_const - 1.0_rt) * state.rho * state.e;
        state.gam1 = gamma_const;

        break;

    case eos_input_rh:

        // dens, enthalpy, and xmass are inputs

#if (!(defined(AMREX_USE_ACC) || defined(AMREX_USE_CUDA)))
        amrex::Error("EOS: eos_input_rh is not supported in this EOS.");
#endif

        break;

    case eos_input_tp:

        // temp, pres, and xmass are inputs

#if (!(defined(AMREX_USE_ACC) || defined(AMREX_USE_CUDA)))
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
        state.cs = std::sqrt(gamma_const * poverrho);

        state.dpdr_e = poverrho;
        state.dpde = (gamma_const - 1.0_rt) * state.rho;

        // Try to avoid the expensive log function.  Since we don't need entropy
        // in hydro solver, set it to an invalid but "nice" value for the plotfile.
        state.s = 1.0_rt;

        break;

    case eos_input_ps:

        // pressure entropy, and xmass are inputs

#if (!(defined(AMREX_USE_ACC) || defined(AMREX_USE_CUDA)))
        amrex::Error("EOS: eos_input_ps is not supported in this EOS.");
#endif

        break;

    case (eos_input_ph):

        // pressure, enthalpy and xmass are inputs

#if (!(defined(AMREX_USE_ACC) || defined(AMREX_USE_CUDA)))
        amrex::Error("EOS: eos_input_ph is not supported in this EOS.");
#endif

        break;

    case eos_input_th:

        // temperature, enthalpy and xmass are inputs

        // This system is underconstrained.

#if (!(defined(AMREX_USE_ACC) || defined(AMREX_USE_CUDA)))
        amrex::Error("EOS: eos_input_th is not a valid input for the gamma law EOS.");
#endif

        break;

    default:

#if (!(defined(AMREX_USE_ACC) || defined(AMREX_USE_CUDA)))
        amrex::Error("EOS: invalid input.");
#endif

        break;

    }

    // Give dpdr a value for the purposes of the composition_derivatives routine.

    state.dpdr = 0.0_rt;

}

#endif
