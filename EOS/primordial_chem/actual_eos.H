#ifndef _actual_eos_H_
#define _actual_eos_H_

// This is a multi-gamma EOS.  Each species can have its own gamma, but
// otherwise, they all act as an ideal gas.

// Note: at the moment, it is not clear what the proper expression for
// a multi-gamma entropy should be, so do not rely on the entropy.

#include <AMReX.H>
#include <network.H>
#include <fundamental_constants.H>
#include <extern_parameters.H>
#include <cmath>
#include <actual_eos_data.H>

const std::string eos_name = "multigamma";

inline
void actual_eos_init ()
{

    // Set the gammas for the species -- we have some runtime parameters
    // that can override the default gammas for a few named species.

    for (int n = 0; n < NumSpec; ++n) {
        gammas[n] = eos_gamma_default;
        spmasses[n] = 1.672621637e-24;
    }

    int idx;

    idx = network_spec_index(species_a_name);
    if (idx >= 0) {
        gammas[idx] = species_a_gamma;
        spmasses[idx] = species_a_mass;
    }

    idx = network_spec_index(species_b_name);
    if (idx >= 0) {
        gammas[idx] = species_b_gamma;
        spmasses[idx] = species_b_mass;
   }

    idx = network_spec_index(species_c_name);
    if (idx >= 0) {
        gammas[idx] = species_c_gamma;
        spmasses[idx] = species_c_mass;
    }

    idx = network_spec_index(species_d_name);
    if (idx >= 0) {
        gammas[idx] = species_d_gamma;
        spmasses[idx] = species_d_mass;
    }

    idx = network_spec_index(species_e_name);
    if (idx >= 0) {
        gammas[idx] = species_e_gamma;
        spmasses[idx] = species_e_mass;
    }

    idx = network_spec_index(species_f_name);
    if (idx >= 0) {
        gammas[idx] = species_f_gamma;
        spmasses[idx] = species_f_mass;
    }

    idx = network_spec_index(species_g_name);
    if (idx >= 0) {
        gammas[idx] = species_g_gamma;
        spmasses[idx] = species_g_mass;
    }

    idx = network_spec_index(species_h_name);
    if (idx >= 0) {
        gammas[idx] = species_h_gamma;
        spmasses[idx] = species_h_mass;
    }

    idx = network_spec_index(species_i_name);
    if (idx >= 0) {
        gammas[idx] = species_i_gamma;
        spmasses[idx] = species_i_mass;
    }

    idx = network_spec_index(species_j_name);
    if (idx >= 0) {
        gammas[idx] = species_j_gamma;
        spmasses[idx] = species_j_mass;
    }

    idx = network_spec_index(species_k_name);
    if (idx >= 0) {
        gammas[idx] = species_k_gamma;
        spmasses[idx] = species_k_mass;
    }

    idx = network_spec_index(species_l_name);
    if (idx >= 0) {
        gammas[idx] = species_l_gamma;
        spmasses[idx] = species_l_mass;
    }

    idx = network_spec_index(species_m_name);
    if (idx >= 0) {
        gammas[idx] = species_m_gamma;
        spmasses[idx] = species_m_mass;
    }

    idx = network_spec_index(species_n_name);
    if (idx >= 0) {
        gammas[idx] = species_n_gamma;
        spmasses[idx] = species_n_mass;
    }

}


template <typename I>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
bool is_input_valid (I input)
{
  static_assert(std::is_same<I, eos_input_t>::value, "input must be an eos_input_t");

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
    static_assert(std::is_same<I, eos_input_t>::value, "input must be an eos_input_t");

    // Get the mass of a nucleon from Avogadro's number.
    const Real m_nucleon = 1.0_rt / C::n_A;
    const Real gasconstant = 8.3144725e7;
    const Real protonmass = 1.672621637e-24;
    // Special gamma factors
    Real sum_Abarinv = 0.0_rt;
    Real sum_gammasinv = 0.0_rt;
    for (int n = 0; n < NumSpec; ++n) {
        sum_Abarinv = sum_Abarinv + (state.xn[n]/(spmasses[n]/protonmass));
        sum_gammasinv = sum_gammasinv + (state.xn[n] / ((spmasses[n]/protonmass)*(gammas[n]-1.0)) );
    }

    sum_gammasinv = sum_gammasinv/sum_Abarinv;

    //-------------------------------------------------------------------------
    // For all EOS input modes EXCEPT eos_input_rt, first compute dens
    // and temp as needed from the inputs.
    //-------------------------------------------------------------------------

    Real temp, dens, eint;

    switch (input) {

    case eos_input_rt:

        // dens, temp and xmass are inputs

        // We don't need to do anything here
        temp = state.T;
        dens = state.rho;
        eint = sum_gammasinv * sum_Abarinv * gasconstant * state.T;
        break;

    case eos_input_re:

        // dens, energy, and xmass are inputs

        // Solve for the temperature
        // e = k T / [(mu m_nucleon)*(gamma-1)]
        if constexpr (has_energy<T>::value) {
            dens = state.rho;
            temp = state.e /( sum_gammasinv * gasconstant * sum_Abarinv);
            eint = state.e;
        }

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

    Real pres = dens * temp * gasconstant * sum_Abarinv;
    Real ener = eint;
    Real enth = ener + pres / dens;

    // Compute the pressure simply from the ideal gas law, and the
    // specific internal energy using the gamma-law EOS relation.
    if constexpr (has_pressure<T>::value) {
        state.p = pres;
    }

    if constexpr (has_energy<T>::value) {
        state.e = eint;
    }

    // Enthalpy is h = e + p/rho
    if constexpr (has_enthalpy<T>::value) {
        state.h = enth;
    }

    // Compute the thermodynamic derivatives and specific heats
    Real dpdT = pres / temp;
    Real dpdr = pres / dens;
    if constexpr (has_pressure<T>::value) {
        state.dpdT = dpdT;
        state.dpdr = dpdr;
    }

    Real dedT = ener / temp;
    Real dedr = 0.0_rt;
    if constexpr (has_energy<T>::value) {
        state.dedT = dedT;
        state.dedr = dedr;
    }

    Real dsdT = 0.0_rt;
    Real dsdr = 0.0_rt;
    if constexpr (has_entropy<T>::value) {
        state.dsdT = dsdT;
        state.dsdr = dsdr;
    }

    Real dhdT = dedT + dpdT / dens;
    Real dhdr = 0.0_rt;
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

#ifdef EXTRA_THERMO
    // These need to be worked out.
    if constexpr (has_pressure<T>::value) {
        state.dpdA = 0.0_rt;
        state.dpdZ = 0.0_rt;
    }

    if constexpr (has_energy<T>::value) {
        state.dedA = 0.0_rt;
        state.dedZ = 0.0_rt;
    }
#endif

}



inline
void actual_eos_finalize ()
{
}

#endif
