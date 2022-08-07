#ifndef _actual_eos_H_
#define _actual_eos_H_

// This is a mult-gamma EOS for primordial ISM chemistry.
// Each species can have its own gamma, but
// otherwise, they all act as an ideal gas.

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

    // Set the gammas & masses for the species

    for (int n = 0; n < NumSpec; ++n) {
        gammas[n] = eos_gamma_default;
        spmasses[n] = 1.67353251819e-24;
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

    const Real gasconstant = 8.3144725e7;
    // Special gamma factors
    Real sum_Abarinv = 0.0_rt;
    Real sum_gammasinv = 0.0_rt;
    Real rhotot = 0.0_rt;

    for (int n = 0; n < NumSpec; ++n) {
        rhotot += state.xn[n]*spmasses[n];
    }

    for (int n = 0; n < NumSpec; ++n) {
<<<<<<< HEAD
        sum_Abarinv = sum_Abarinv + (state.xn[n]/spmasses[n]);
        sum_gammasinv = sum_gammasinv + (state.xn[n] / (spmasses[n]*(gammas[n]-1.0)) );
    }

=======
        sum_Abarinv += state.xn[n];
        sum_gammasinv += (state.xn[n]*protonmass/rhotot) * (1.0/(gammas[n]-1.0));
    }

    sum_Abarinv *= protonmass/rhotot;
    sum_gammasinv /= sum_Abarinv;
>>>>>>> f9a80208 (cleaned up EOS)

    //-------------------------------------------------------------------------
    // For all EOS input modes EXCEPT eos_input_rt, first compute dens
    // and temp as needed from the inputs.
    //-------------------------------------------------------------------------

    Real temp, dens;

    switch (input) {

    case eos_input_rt:

        // dens, temp and xmass are inputs
        // We don't need to do anything here

        temp = state.T;
        dens = state.rho;

        break;

    case eos_input_re:

        // dens, energy, and xmass are inputs
        // Solve for the temperature

        if constexpr (has_energy<T>::value) {
            dens = state.rho;
<<<<<<< HEAD
            temp = state.e * sum_gammasinv * gasconstant / sum_Abarinv;
=======
            // stop the integration if the internal energy < 0
            AMREX_ASSERT(state.e > 0.);
            temp = state.e /( sum_gammasinv * gasconstant * sum_Abarinv);
            eint = state.e;
>>>>>>> f9a80208 (cleaned up EOS)
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

<<<<<<< HEAD
    Real pres = dens * temp * gasconstant * sum_Abarinv;
    Real ener = 0.0_rt; //currently undefined, as we need kinetic energy total define total energy
    Real enth = ener + pres / dens;

    // Compute the pressure simply from the ideal gas law, and the
    // specific internal energy using the gamma-law EOS relation.
    if constexpr (has_pressure<T>::value) {
        state.p = pres;
    }
=======
>>>>>>> f9a80208 (cleaned up EOS)
    if constexpr (has_energy<T>::value) {
        state.e = ener;
    }

    Real dedT = sum_gammasinv * sum_Abarinv * gasconstant;
    if constexpr (has_energy<T>::value) {
        state.dedT = dedT;
    }

    //if constexpr (has_energy<T>::value) {
    //    state.cv = state.dedT;
    //}

}



inline
void actual_eos_finalize ()
{
}

#endif
