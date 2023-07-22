#ifndef ACTUAL_NETWORK_H
#define ACTUAL_NETWORK_H

#define NEW_NETWORK_IMPLEMENTATION

#include <AMReX_REAL.H>

using namespace amrex;

#include <network_properties.H>
#include <extern_parameters.H>
#include <rhs_type.H>

AMREX_INLINE
void actual_network_init() {}

const std::string network_name = "powerlaw";

namespace Rates {
    enum NetworkRates {
        powerlaw_rate = 1,
        NumRates = powerlaw_rate
    };
}

namespace RHS {

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    constexpr rhs_t rhs_data (int rate)
    {
        using namespace Species;
        using namespace Rates;

        rhs_t data;

        switch (rate) {

        case powerlaw_rate:
            data.species_A = Fuel;
            data.species_D = Ash;

            data.number_A = 2;
            data.number_D = 1;

            data.exponent_A = 2;
            data.exponent_D = 1;

            data.apply_identical_particle_factor = 0;
            break;

        }

        return data;
    }

    template<int rate>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void evaluate_analytical_rate (const rhs_state_t& state, rate_t& rates)
    {
        using namespace Species;
        using namespace Rates;

        if constexpr (rate == powerlaw_rate) {
            if (state.tf.temp < network_rp::f_act * network_rp::T_burn_ref) {
                rates.fr = 0.0_rt;
            }
            else {
                rates.fr = network_rp::rtilde / network_rp::rho_burn_ref *
                           std::pow(state.tf.temp / network_rp::T_burn_ref, network_rp::nu);
            }

            rates.rr = 0.0_rt;
        }
    }

    template<int rate>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void postprocess_rate (const rhs_state_t& state, rate_t& rates,
                           rate_t& rates1, rate_t& rates2, rate_t& rates3)
    {
        // Nothing to do for this network.
    }

    template<int spec>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real ener_gener_rate (const rhs_state_t& rhs_state, Real const& dydt)
    {
        if constexpr (spec == 2) {
            return dydt * NetworkProperties::aion(spec) * network_rp::specific_q_burn;
        }
        else {
            return 0.0;
        }
    }

} // namespace RHS

#endif
