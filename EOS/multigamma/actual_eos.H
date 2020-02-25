#ifndef _actual_eos_H_
#define _actual_eos_H_

// This is a multi-gamma EOS.  Each species can have its own gamma, but
// otherwise, they all act as an ideal gas.

// Note: at the moment, it is not clear what the proper expression for
// a multi-gamma entropy should be, so do not rely on the entropy.

#include <AMReX.H>
#include <network.H>
#include <actual_eos_data.H>
#include <fundamental_constants.H>
#include <extern_parameters.H>
#include <cmath>

const std::string eos_name = "multigamma";

inline
void actual_eos_init ()
{

    // Set the gammas for the species -- we have some runtime parameters
    // that can override the default gammas for a few named species.

    for (int n = 0; n < NumSpec; ++n) {
        gammas[n] = eos_gamma_default;
    }

    int idx;

    idx = network_spec_index(species_a_name);
    if (idx >= 0) {
        gammas[idx] = species_a_gamma;
    }

    idx = network_spec_index(species_b_name);
    if (idx >= 0) {
        gammas[idx] = species_b_gamma;
    }

    idx = network_spec_index(species_c_name);
    if (idx >= 0) {
        gammas[idx] = species_c_gamma;
    }

}



AMREX_GPU_HOST_DEVICE inline
void actual_eos (eos_input_t input, eos_t& state)
{
    // Get the mass of a nucleon from Avogadro's number.
    const Real m_nucleon = 1.0_rt / n_A;

    // Special gamma factors
    Real sumY_gm1 = 0.0_rt;
    Real sumYg_gm1 = 0.0_rt;
    for (int n = 0; n < NumSpec; ++n) {
        sumY_gm1 = sumY_gm1 + state.xn[n] * aion_inv[n] / (gammas[n] - 1.0_rt);
        sumYg_gm1 = sumYg_gm1 + state.xn[n] * gammas[n] * aion_inv[n] / (gammas[n] - 1.0_rt);
    }

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

    case eos_input_rh:

        // dens, enthalpy, and xmass are inputs

        // Solve for the temperature:
        // h = e + p/rho = (p/rho)*[1 + 1/(gamma-1)] = (p/rho)*gamma/(gamma-1)
        dens = state.rho;
        temp = state.h * (m_nucleon / k_B) / sumYg_gm1;

        break;

    case eos_input_tp:

        // temp, pres, and xmass are inputs

        // Solve for the density:
        // p = rho k T / (abar m_nucleon)
        dens = state.p * state.abar * (m_nucleon / k_B) / state.T;
        temp = state.T;

        break;

    case eos_input_rp:

        // dens, pres, and xmass are inputs

        // Solve for the temperature:
        // p = rho k T / (mu m_nucleon)
        dens = state.rho;
        temp = state.p * state.abar * (m_nucleon / k_B) / state.rho;

        break;

    case eos_input_re:

        // dens, energy, and xmass are inputs

        // Solve for the temperature
        // e = k T / [(mu m_nucleon)*(gamma-1)]
        dens = state.rho;
        temp = state.e * (m_nucleon / k_B) / sumY_gm1;

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
        dens = state.p * state.abar / state.h * sumYg_gm1;
        temp = state.p * state.abar * (m_nucleon / k_B) / dens;

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

    // Compute the pressure simply from the ideal gas law, and the
    // specific internal energy using the gamma-law EOS relation.
    state.p = dens * (k_B / m_nucleon) * temp / state.abar;
    state.e = (k_B / m_nucleon) * temp * sumY_gm1;

    // Enthalpy is h = e + p/rho
    state.h = state.e + state.p / dens;

    // entropy (per gram) -- this is wrong. Not sure what the expression
    // is for a multigamma gas
    state.s = ((k_B / m_nucleon) / state.abar) *
              (2.5_rt + std::log((std::pow(state.abar * m_nucleon, 2.5_rt) / dens) *
                                 std::pow(k_B * temp, 1.5_rt) /
                                 std::pow(2.0_rt * M_PI * hbar * hbar, 1.5_rt)));

    // Compute the thermodynamic derivatives and specific heats
    state.dpdT = state.p / temp;
    state.dpdr = state.p / dens;
    state.dedT = state.e / temp;
    state.dedr = 0.0_rt;
    state.dsdT = 0.0_rt;
    state.dsdr = 0.0_rt;
    state.dhdT = state.dedT + state.dpdT / dens;
    state.dhdr = 0.0_rt;

    state.cv = state.dedT;
    state.cp = state.h / state.T;

    state.gam1 = state.cp / state.cv;

    state.dpdr_e = state.dpdr - state.dpdT * state.dedr / state.dedT;
    state.dpde   = state.dpdT / state.dedT;

#ifdef EXTRA_THERMO
    // These need to be worked out.
    state.dpdA = 0.0_rt;
    state.dedA = 0.0_rt;

    state.dpdZ = 0.0_rt;
    state.dedZ = 0.0_rt;
#endif

    // Sound speed
    state.cs = std::sqrt(state.gam1 * state.p / dens);

}



inline
void actual_eos_finalize ()
{
}

#endif
