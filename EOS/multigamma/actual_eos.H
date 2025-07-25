#ifndef ACTUAL_EOS_H
#define ACTUAL_EOS_H

// This is a multi-gamma EOS.  Each species can have its own gamma, but
// otherwise, they all act as an ideal gas.

// Note: at the moment, it is not clear what the proper expression for
// a multi-gamma entropy should be, so do not rely on the entropy.

#include <string_view>
#include <AMReX.H>
#include <network.H>
#include <fundamental_constants.H>
#include <extern_parameters.H>
#include <cmath>
#include <actual_eos_data.H>

constexpr std::string_view eos_name = "multigamma";

inline
void actual_eos_init ()
{

    // Set the gammas for the species -- we have some runtime parameters
    // that can override the default gammas for a few named species.

    for (int n = 0; n < NumSpec; ++n) {
        gammas[n] = eos_rp::eos_gamma_default;
    }

    int idx;

    idx = network_spec_index(eos_rp::species_a_name);
    if (idx >= 0) {
        gammas[idx] = eos_rp::species_a_gamma;
    }

    idx = network_spec_index(eos_rp::species_b_name);
    if (idx >= 0) {
        gammas[idx] = eos_rp::species_b_gamma;
    }

    idx = network_spec_index(eos_rp::species_c_name);
    if (idx >= 0) {
        gammas[idx] = eos_rp::species_c_gamma;
    }

}


template <typename I>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
bool is_input_valid (I input)
{
  static_assert(std::is_same_v<I, eos_input_t>, "input must be an eos_input_t");

  bool valid = true;

  if (input == eos_input_ps ||
      input == eos_input_th) {
    valid = false;
  }

  return valid;
}


template <typename I, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_eos (I input, T& state)
{
    static_assert(std::is_same_v<I, eos_input_t>, "input must be an eos_input_t");

    // Get the mass of a nucleon from Avogadro's number.
    const amrex::Real m_nucleon = 1.0_rt / C::n_A;

    // Special gamma factors
    amrex::Real sumY_gm1 = 0.0_rt;
    amrex::Real sumYg_gm1 = 0.0_rt;
    for (int n = 0; n < NumSpec; ++n) {
        sumY_gm1 = sumY_gm1 + state.xn[n] * aion_inv[n] / (gammas[n] - 1.0_rt);
        sumYg_gm1 = sumYg_gm1 + state.xn[n] * gammas[n] * aion_inv[n] / (gammas[n] - 1.0_rt);
    }

    //-------------------------------------------------------------------------
    // For all EOS input modes EXCEPT eos_input_rt, first compute dens
    // and temp as needed from the inputs.
    //-------------------------------------------------------------------------

    amrex::Real temp, dens;

    switch (input) {

    case eos_input_rt:

        // dens, temp and xmass are inputs

        // We don't need to do anything here
        temp = state.T;
        dens = state.rho;

        break;

    case eos_input_rh:

        // dens, enthalpy, and xmass are inputs

        // Solve for the temperature:
        // h = e + p/rho = (p/rho)*[1 + 1/(gamma-1)] = (p/rho)*gamma/(gamma-1)
        if constexpr (has_enthalpy<T>::value) {
            dens = state.rho;
            temp = state.h * (m_nucleon / C::k_B) / sumYg_gm1;
        }

        break;

    case eos_input_tp:

        // temp, pres, and xmass are inputs

        // Solve for the density:
        // p = rho k T / (abar m_nucleon)
        if constexpr (has_pressure<T>::value) {
            dens = state.p * state.abar * (m_nucleon / C::k_B) / state.T;
            temp = state.T;
        }

        break;

    case eos_input_rp:

        // dens, pres, and xmass are inputs

        // Solve for the temperature:
        // p = rho k T / (mu m_nucleon)
        if constexpr (has_pressure<T>::value) {
            dens = state.rho;
            temp = state.p * state.abar * (m_nucleon / C::k_B) / state.rho;
        }

        break;

    case eos_input_re:

        // dens, energy, and xmass are inputs

        // Solve for the temperature
        // e = k T / [(mu m_nucleon)*(gamma-1)]
        if constexpr (has_energy<T>::value) {
            dens = state.rho;
            temp = state.e * (m_nucleon / C::k_B) / sumY_gm1;
        }

        break;

    case eos_input_ps:

        // pressure entropy, and xmass are inputs
#ifndef AMREX_USE_GPU
        amrex::Error("eos_input_ps is not supported");
#endif

        break;

    case eos_input_ph:

        // pressure, enthalpy and xmass are inputs

        // Solve for temperature and density
        if constexpr (has_pressure<T>::value && has_enthalpy<T>::value) {
            dens = state.p * state.abar / state.h * sumYg_gm1;
            temp = state.p * state.abar * (m_nucleon / C::k_B) / dens;
        }

        break;

    case eos_input_th:

        // temperature, enthalpy and xmass are inputs
#ifndef AMREX_USE_GPU
        amrex::Error("eos_input_th is not supported");
#endif

        break;

    default:

#ifndef AMREX_USE_GPU
        amrex::Error("EOS: invalid input.");
#endif

        break;

    }

    //-------------------------------------------------------------------------
    // Now we have the density and temperature (and mass fractions /
    // mu), regardless of the inputs.
    //-------------------------------------------------------------------------

    state.T   = temp;
    state.rho = dens;
    amrex::Real pres = dens * (C::k_B / m_nucleon) * temp / state.abar;
    amrex::Real ener = (C::k_B / m_nucleon) * temp * sumY_gm1;
    amrex::Real enth = ener + pres / dens;

    // Compute the pressure simply from the ideal gas law, and the
    // specific internal energy using the gamma-law EOS relation.
    if constexpr (has_pressure<T>::value) {
        state.p = pres;
    }
    if constexpr (has_energy<T>::value) {
        state.e = ener;
    }

    // Enthalpy is h = e + p/rho
    if constexpr (has_enthalpy<T>::value) {
        state.h = enth;
    }

    // entropy (per gram) -- this is wrong. Not sure what the expression
    // is for a multigamma gas
    if constexpr (has_entropy<T>::value) {
        state.s = ((C::k_B / m_nucleon) / state.abar) *
                  (2.5_rt + std::log((std::pow(state.abar * m_nucleon, 2.5_rt) / dens) *
                                     std::pow(C::k_B * temp, 1.5_rt) /
                                     std::pow(2.0_rt * M_PI * C::hbar * C::hbar, 1.5_rt)));
    }

    // Compute the thermodynamic derivatives and specific heats
    amrex::Real dpdT = pres / temp;
    amrex::Real dpdr = pres / dens;
    if constexpr (has_pressure<T>::value) {
        state.dpdT = dpdT;
        state.dpdr = dpdr;
    }

    amrex::Real dedT = ener / temp;
    amrex::Real dedr = 0.0_rt;
    if constexpr (has_energy<T>::value) {
        state.dedT = dedT;
        state.dedr = dedr;
    }

    amrex::Real dsdT = 0.0_rt;
    amrex::Real dsdr = 0.0_rt;
    if constexpr (has_entropy<T>::value) {
        state.dsdT = dsdT;
        state.dsdr = dsdr;
    }

    amrex::Real dhdT = dedT + dpdT / dens;
    amrex::Real dhdr = 0.0_rt;
    if constexpr (has_enthalpy<T>::value) {
        state.dhdT = dhdT;
        state.dhdr = dhdr;
    }

    if constexpr (has_energy<T>::value) {
        state.cv = state.dedT;
    }
    if constexpr (has_pressure<T>::value && has_energy<T>::value) {
        state.cp = enth / state.T;

        state.gam1 = state.cp / state.cv;

        state.dpdr_e = state.dpdr - state.dpdT * state.dedr / state.dedT;
        state.dpde   = state.dpdT / state.dedT;
    }

    // These need to be worked out.
    if constexpr (has_dpdA<T>::value) {
        state.dpdA = 0.0_rt;
    }
    if constexpr (has_dpdZ<T>::value) {
        state.dpdZ = 0.0_rt;
    }

    if constexpr (has_dedA<T>::value) {
        state.dedA = 0.0_rt;
    }
    if constexpr (has_dedZ<T>::value) {
        state.dedZ = 0.0_rt;
    }

    // Sound speed
    if constexpr (has_pressure<T>::value && has_energy<T>::value) {
        state.cs = std::sqrt(state.gam1 * state.p / dens);
    }
}



inline
void actual_eos_finalize ()
{
}

#endif
