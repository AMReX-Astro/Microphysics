#ifndef ACTUAL_NETWORK_H
#define ACTUAL_NETWORK_H

#define NEW_NETWORK_IMPLEMENTATION

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>

using namespace amrex;

#include <fundamental_constants.H>
#include <network_properties.H>
#include <rhs_type.H>
#include <extern_parameters.H>
#include <tfactors.H>
#include <aprox_rates.H>

void actual_network_init();

const std::string network_name = "iso7";

namespace network
{
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> bion;
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> mion;
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> wion;
}

namespace Rates {
    enum NetworkRates {
        C12_He4_to_O16 = 1,
        He4_He4_He4_to_C12,
        C12_C12_to_Ne20_He4,
        C12_O16_to_Mg24_He4,
        C12_O16_to_Si28,
        O16_O16_to_Si28_He4,
        O16_He4_to_Ne20,
        Ne20_He4_to_Mg24,
        Mg24_He4_to_Si28,
        Ca40_He4_to_Ti44,
        Si28_7He4_to_Ni56,
        NumRates = Si28_7He4_to_Ni56
    };
}

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

        data.forward_branching_ratio = 1.0_rt;
        data.reverse_branching_ratio = 1.0_rt;

        data.apply_identical_particle_factor = 1;

        data.rate_can_be_tabulated = 1;

        data.screen_forward_reaction = 1;
        data.screen_reverse_reaction = 0; // For consistency with original iso7

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
            break;

        case C12_C12_to_Ne20_He4:
            data.species_A = C12;
            data.species_D = Ne20;
            data.species_E = He4;

            data.number_A = 2;
            data.number_D = 1;
            data.number_E = 1;

            data.exponent_A = 2;
            data.exponent_D = 1;
            data.exponent_E = 1;
            break;

        case C12_O16_to_Mg24_He4:
            data.species_A = C12;
            data.species_B = O16;
            data.species_D = Mg24;
            data.species_E = He4;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;
            data.number_E = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            data.exponent_E = 1;

            // Equal probability branching ratio between (Mg24 + He4) and (Si28) endpoints
            data.forward_branching_ratio = 0.5_rt;
            break;

        case C12_O16_to_Si28:
            data.species_A = C12;
            data.species_B = O16;
            data.species_D = Si28;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;

            // Equal probability branching ratio between (Mg24 + He4) and (Si28) endpoints
            data.forward_branching_ratio = 0.5_rt;
            break;

        case O16_O16_to_Si28_He4:
            data.species_A = O16;
            data.species_D = Si28;
            data.species_E = He4;

            data.number_A = 2;
            data.number_D = 1;
            data.number_E = 1;

            data.exponent_A = 2;
            data.exponent_D = 1;
            data.exponent_E = 1;
            break;

        case O16_He4_to_Ne20:
            data.species_A = O16;
            data.species_B = He4;
            data.species_D = Ne20;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            break;

        case Ne20_He4_to_Mg24:
            data.species_A = Ne20;
            data.species_B = He4;
            data.species_D = Mg24;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            break;

        case Mg24_He4_to_Si28:
            data.species_A = Mg24;
            data.species_B = He4;
            data.species_D = Si28;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            break;

        case Ca40_He4_to_Ti44:
            data.species_A = Ca40;
            data.species_B = He4;
            data.species_D = Ti44;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            break;

        case Si28_7He4_to_Ni56:
            data.species_A = Si28;
            data.species_B = He4;
            data.species_D = Ni56;

            data.number_A = 1;
            data.number_B = 7;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;

            data.rate_can_be_tabulated = 0;

            data.additional_reaction_1 = Ca40_He4_to_Ti44;
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

        if constexpr (rate == C12_He4_to_O16) {
            rate_c12ag(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == He4_He4_He4_to_C12) {
            rate_triplealf(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == C12_C12_to_Ne20_He4) {
            rate_c12c12(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == C12_O16_to_Mg24_He4) {
            rate_c12o16(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == C12_O16_to_Si28) {
            rate_c12o16(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == O16_O16_to_Si28_He4) {
            rate_o16o16(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == O16_He4_to_Ne20) {
            rate_o16ag(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == Ne20_He4_to_Mg24) {
            rate_ne20ag(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == Mg24_He4_to_Si28) {
            rate_mg24ag(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == Ca40_He4_to_Ti44) {
            rate_ca40ag(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
    }

    template<int rate>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void postprocess_rate (const rhs_state_t& state, rate_t& rates,
                           rate_t& rates1, [[maybe_unused]] rate_t& rates2, [[maybe_unused]] rate_t& rates3)
    {
        using namespace Species;
        using namespace Rates;

        if constexpr (rate == Si28_7He4_to_Ni56) {
            // first rate corresponds to Ca40_He4_to_Ti44

            if (state.tf.t9 > 2.5_rt && state.y(C12) + state.y(O16) <= 4.0e-3_rt) {

                Real t992  = state.tf.t972 * state.tf.t9;
                Real t9i92 = 1.0_rt / t992;

                Real yeff_ca40   = t9i92 * std::exp(239.42_rt * state.tf.t9i - 74.741_rt);
                Real yeff_ca40dt = -yeff_ca40 * (239.42_rt * state.tf.t9i2 + 4.5_rt * state.tf.t9i);

                Real yeff_ti44   = t992  * std::exp(-274.12_rt * state.tf.t9i + 74.914_rt);
                Real yeff_ti44dt = yeff_ti44*(274.12_rt * state.tf.t9i2 + 4.5_rt * state.tf.t9i);

                Real denom = std::pow(state.rho * state.y(He4), 3.0_rt);

                rates.fr   = yeff_ca40 * denom * rates1.fr;
                rates.frdt = (yeff_ca40dt * rates1.fr + yeff_ca40 * rates1.frdt) * denom * 1.0e-9_rt;

                Real zz = 1.0_rt / denom;
                rates.rr = amrex::min(1.0e10_rt, yeff_ti44 * rates1.rr * zz);

                if (rates.rr == 1.0e10_rt) {
                    rates.rrdt = 0.0_rt;
                } else {
                    rates.rrdt = (yeff_ti44dt * rates1.rr + yeff_ti44 * rates1.rrdt) * zz * 1.0e-9_rt;
                }
            }
        }
    }

    template<int spec>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real ener_gener_rate (Real const& dydt)
    {
        return dydt * network::mion(spec) * C::Legacy::enuc_conv2;
    }

} // namespace RHS

#endif
