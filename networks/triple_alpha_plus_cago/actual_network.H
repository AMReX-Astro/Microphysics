#ifndef actual_network_H
#define actual_network_H

#define NEW_NETWORK_IMPLEMENTATION

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>

#include <fundamental_constants.H>
#include <network_properties.H>
#include <rhs_type.H>
#include <aprox_rates.H>

using namespace amrex;

AMREX_INLINE
void actual_network_init() {}

const std::string network_name = "triple_alpha_plus_cago";

namespace network
{
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real bion (int spec) {
        using namespace Species;

        // our convention is that binding energy is negative.  The following are
        // the binding energies per unit mass (erg / g) obtained by converting
        // the energies in MeV to erg then multiplying by (N_A / aion) where
        // N_A = 6.0221415e23 is Avogadro's number

        if (spec == He4) {
            return 28.29603_rt;
        }
        else if (spec == C12) {
            return 92.16294_rt;
        }
        else if (spec == O16) {
            return 127.62093_rt;
        }
        else if (spec == Fe56) {
            return 492.25389_rt;
        }

        // Return zero if we don't recognize the species.

        return 0.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real mion (int spec) {
        Real A = NetworkProperties::aion(spec);
        Real Z = NetworkProperties::zion(spec);

        return (A - Z) * C::Legacy::m_n + Z * (C::Legacy::m_p + C::Legacy::m_e) - bion(spec) * C::Legacy::MeV2gr;
    }
}

namespace Rates
{
    enum NetworkRates {
        He4_He4_He4_to_C12 = 1,
        C12_He4_to_O16,
        NumRates = C12_He4_to_O16
    };
};

namespace RHS {
    
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    constexpr rhs_t rhs_data (int rate)
    {
        using namespace Species;
        using namespace Rates;

        rhs_t data{};

        data.species_A = -1;
        data.species_B = -1;
        data.species_C = -1;
        data.species_D = -1;
        data.species_E = -1;
        data.species_F = -1;

        data.number_A = 0;
        data.number_B = 0;
        data.number_C = 0;
        data.number_D = 0;
        data.number_E = 0;
        data.number_F = 0;

        data.exponent_A = 0;
        data.exponent_B = 0;
        data.exponent_C = 0;
        data.exponent_D = 0;
        data.exponent_E = 0;
        data.exponent_F = 0;

        data.forward_branching_ratio = 1.0;
        data.reverse_branching_ratio = 1.0;

        data.apply_identical_particle_factor = 1;

        data.rate_can_be_tabulated = 1;

        data.screen_forward_reaction = 1;
        data.screen_reverse_reaction = 1;

        data.additional_reaction_1 = -1;
        data.additional_reaction_2 = -1;
        data.additional_reaction_3 = -1;

        switch (rate) {

        case He4_He4_He4_to_C12:
            data.species_A = He4;
            data.species_D = C12;

            data.number_A = 3;
            data.number_D = 1;

            data.exponent_A = 3;
            data.exponent_D = 1;

            // This network doesn't allow the reverse reaction
            data.reverse_branching_ratio = 0.0_rt;
            break;

        case C12_He4_to_O16:
            data.species_A = C12;
            data.species_B = He4;
            data.species_D = O16;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;

            // This network doesn't allow the reverse reaction
            data.reverse_branching_ratio = 0.0_rt;
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

        if constexpr (rate == He4_He4_He4_to_C12) {
            rate_triplealf(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == C12_He4_to_O16) {
            rate_c12ag(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
    }

    template<int rate>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void postprocess_rate ([[maybe_unused]] const rhs_state_t& state, rate_t& rates,
                           rate_t& rates1, rate_t& rates2, rate_t& rates3)
    {
        using namespace Species;
        using namespace Rates;

        // Nothing to do for this network.
    }

    template<int spec>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real ener_gener_rate (Real const& dydt)
    {
        return dydt * network::mion(spec) * C::Legacy::enuc_conv2;
    }

} // namespace RHS

#endif
