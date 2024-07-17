#ifndef actual_network_H
#define actual_network_H

#define NEW_NETWORK_IMPLEMENTATION

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>

#include <fundamental_constants.H>
#include <network_properties.H>
#include <rhs_type.H>

using namespace amrex::literals;

AMREX_INLINE
void actual_network_init() {}

const std::string network_name = "ignition_chamulak";

namespace Rates
{
    enum NetworkRates {
        chamulak_rate = 1,
        NumRates = chamulak_rate
    };
};

namespace RHS {

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    constexpr rhs_t rhs_data (int rate)
    {
        using namespace Species;
        using namespace Rates;

        rhs_t data;

        switch (rate) {

        case chamulak_rate:
            data.species_A = C12;
            data.species_D = Ash;

            // The composite "ash" particle we are creating has a A = 18,
            // while we are consuming 2 C12 nucleons to create it. In order
            // to conserve mass, we need to scale the number of ash particles
            // created by 24 / 18 = 4 / 3.

            data.number_A = 2.0_rt;
            data.number_D = 4.0_rt / 3.0_rt;

            data.exponent_A = 2;
            data.exponent_D = 1;

            // For each reaction, we lose 2.93 C12 nuclei (for a single rate,
            // C12+C12, we would lose 2. In Chamulak et al. (2008), they say a
            // value of 2.93 captures the energetics from a larger network.

            data.forward_branching_ratio = 2.93e0_rt / 2.0_rt;
            break;

        }

        return data;
    }

    template<int rate, typename number_t>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void evaluate_analytical_rate (const rhs_state_t<number_t>& state, rate_t& rates)
    {
        using namespace Species;
        using namespace Rates;

        if constexpr (rate == chamulak_rate) {
            // compute some often used temperature constants
            amrex::Real T9     = state.tf.t9;
            amrex::Real dT9dt  = 1.0_rt / 1.0e9_rt;
            amrex::Real T9a    = T9 / (1.0_rt + 0.0396e0_rt * T9);
            amrex::Real dT9adt = (T9a / T9 - (T9a / (1.0_rt + 0.0396e0_rt * T9)) * 0.0396e0_rt) * dT9dt;

            // compute the CF88 rate
            amrex::Real scratch    = std::pow(T9a, 1.0_rt / 3.0_rt);
            amrex::Real dscratchdt = (1.0_rt / 3.0_rt) * std::pow(T9a, -2.0_rt / 3.0_rt) * dT9adt;

            amrex::Real a    = 4.27e26_rt * std::pow(T9a, 5.0_rt / 6.0_rt) * std::pow(T9, -1.5e0_rt);
            amrex::Real dadt = (5.0_rt / 6.0_rt) * (a / T9a) * dT9adt - 1.5e0_rt * (a / T9) * dT9dt;

            amrex::Real b    = std::exp(-84.165e0_rt / scratch - 2.12e-3_rt * T9 * T9 * T9);
            amrex::Real dbdt = (84.165e0_rt * dscratchdt / (scratch * scratch) - 3.0_rt * 2.12e-3_rt * T9 * T9 * dT9dt) * b;

            rates.fr = a * b;
            rates.frdt = dadt * b + a * dbdt;

            rates.rr = 0.0_rt;
            rates.rrdt = 0.0_rt;
        }
    }

    template<int rate, typename number_t>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void postprocess_rate (const rhs_state_t<number_t>& state, rate_t& rates,
                           rate_t& rates1, rate_t& rates2, rate_t& rates3)
    {
        // Nothing to do for this network.
    }

    template<int spec, typename number_t>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    amrex::Real ener_gener_rate (const rhs_state_t<number_t>& rhs_state, amrex::Real const& dydt)
    {
        // Chamulak et al. provide the q-value resulting from C12 burning,
        // given as 3 different values (corresponding to 3 different densities).
        // Here we do a simple quadratic fit to the 3 values provided (see
        // Chamulak et al., p. 164, column 2).

        // our convention is that the binding energies are negative.  We convert
        // from the MeV values that are traditionally written in astrophysics
        // papers by multiplying by 1.e6 eV/MeV * 1.60217646e-12 erg/eV.  The
        // MeV values are per nucleus, so we divide by aion to make it per
        // nucleon and we multiple by Avogardo's # (6.0221415e23) to get the
        // value in erg/g
        amrex::Real rho9 = rhs_state.rho / 1.0e9_rt;

        // q_eff is effective heat evolved per reaction (given in MeV)
        amrex::Real q_eff = 0.06e0_rt * rho9 * rho9 + 0.02e0_rt * rho9 + 8.83e0_rt;

        // convert from MeV to ergs / gram.  Here 2.93 is the
        // number of C12 nuclei destroyed in a single reaction and 12.0 is
        // the mass of a single C12 nuclei.  Also note that our convention
        // is that binding energies are negative.
        amrex::Real ebin_C12 = -q_eff * C::MeV2eV * C::ev2erg * C::n_A / (2.93e0_rt * 12.0e0_rt);

        if constexpr (spec == Species::C12) {
            return dydt * NetworkProperties::aion(spec) * ebin_C12;
        }
        else {
            return 0.0_rt;
        }
    }

} // namespace RHS

#endif
