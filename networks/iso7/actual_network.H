#ifndef ACTUAL_NETWORK_H
#define ACTUAL_NETWORK_H

#define NEW_NETWORK_IMPLEMENTATION

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>
#include <string_view>

using namespace amrex::literals;

#include <fundamental_constants.H>
#include <network_properties.H>
#include <rhs_type.H>
#include <extern_parameters.H>
#include <tfactors.H>
#include <aprox_rates.H>

AMREX_INLINE
void actual_network_init() {}

constexpr std::string_view network_name = "iso7";

namespace network
{
    template<int spec>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    constexpr amrex::Real bion () {
        using namespace Species;

        static_assert(spec >= 1 && spec <= NumSpec);

        // Set the binding energy of the element
        if constexpr (spec == He4) {
            return 28.29603e0_rt;
        }
        else if constexpr (spec == C12) {
            return 92.16294e0_rt;
        }
        else if constexpr (spec == O16) {
            return 127.62093e0_rt;
        }
        else if constexpr (spec == Ne20) {
            return 160.64788e0_rt;
        }
        else if constexpr (spec == Mg24) {
            return 198.25790e0_rt;
        }
        else if constexpr (spec == Si28) {
            return 236.53790e0_rt;
        }
        else if constexpr (spec == Ni56) {
            return 484.00300e0_rt;
        }

        // Return zero if we don't recognize the species.
        return 0.0_rt;
    }

    template<int spec>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    constexpr amrex::Real mion () {
        static_assert(spec >= 1 && spec <= NumSpec);

        constexpr amrex::Real A = NetworkProperties::aion(spec);
        constexpr amrex::Real Z = NetworkProperties::zion(spec);

        return (A - Z) * C::Legacy::m_n + Z * (C::Legacy::m_p + C::Legacy::m_e) - bion<spec>() * C::Legacy::MeV2gr;
    }

    // Legacy (non-templated) interfaces

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    amrex::Real bion (int spec) {
        using namespace Species;

        amrex::Real b = 0.0_rt;

        // Set the binding energy of the element
        amrex::constexpr_for<1, NumSpec+1>([&] (auto n) {
            if (n == spec) {
                b = bion<n>();
            }
        });

        return b;
    }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    amrex::Real mion (int spec) {
        using namespace Species;

        amrex::Real m = 0.0_rt;

        amrex::constexpr_for<1, NumSpec+1>([&] (auto n) {
            if (n == spec) {
                m = mion<n>();
            }
        });

        return m;
    }
}

#ifdef NSE_NET
namespace NSE_INDEX
{
    // P_index is for photoionization proton
    constexpr int P_index = -1;
    constexpr int H1_index = -1;
    constexpr int N_index = -1;
    constexpr int He4_index = 0;
}
#endif

namespace Rates {
    enum NetworkRates : std::uint8_t {
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

        rhs_t data;

        data.screen_reverse_reaction = 0; // For consistency with original iso7

        switch (rate) {

        case He4_He4_He4_to_C12:
            data.species_A = He4;
            data.species_D = C12;

            data.number_A = 3.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 3;
            data.exponent_D = 1;
            break;

        case C12_He4_to_O16:
            data.species_A = C12;
            data.species_B = He4;
            data.species_D = O16;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            break;

        case C12_C12_to_Ne20_He4:
            data.species_A = C12;
            data.species_D = Ne20;
            data.species_E = He4;

            data.number_A = 2.0_rt;
            data.number_D = 1.0_rt;
            data.number_E = 1.0_rt;

            data.exponent_A = 2;
            data.exponent_D = 1;
            data.exponent_E = 1;
            break;

        case C12_O16_to_Mg24_He4:
            data.species_A = C12;
            data.species_B = O16;
            data.species_D = Mg24;
            data.species_E = He4;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;
            data.number_E = 1.0_rt;

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

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;

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

            data.number_A = 2.0_rt;
            data.number_D = 1.0_rt;
            data.number_E = 1.0_rt;

            data.exponent_A = 2;
            data.exponent_D = 1;
            data.exponent_E = 1;
            break;

        case O16_He4_to_Ne20:
            data.species_A = O16;
            data.species_B = He4;
            data.species_D = Ne20;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            break;

        case Ne20_He4_to_Mg24:
            data.species_A = Ne20;
            data.species_B = He4;
            data.species_D = Mg24;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            break;

        case Mg24_He4_to_Si28:
            data.species_A = Mg24;
            data.species_B = He4;
            data.species_D = Si28;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            break;

        case Ca40_He4_to_Ti44:
            data.species_A = Ca40;
            data.species_B = He4;
            data.species_D = Ti44;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            break;

        case Si28_7He4_to_Ni56:
            data.species_A = Si28;
            data.species_B = He4;
            data.species_D = Ni56;

            data.number_A = 1.0_rt;
            data.number_B = 7.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;

            data.rate_can_be_tabulated = 0;

            data.additional_reaction_1 = Ca40_He4_to_Ti44;
            break;

        }

        return data;
    }

    template<int rate, typename number_t>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void evaluate_analytical_rate (const rhs_state_t<number_t>& state, rate_t<number_t>& rates)
    {
        using namespace Species;
        using namespace Rates;

        if constexpr (rate == C12_He4_to_O16) {
            rate_c12ag(state.tf, 1.0_rt, rates.fr, rates.rr);
        }
        else if constexpr (rate == He4_He4_He4_to_C12) {
            rate_triplealf(state.tf, 1.0_rt, rates.fr, rates.rr);
        }
        else if constexpr (rate == C12_C12_to_Ne20_He4) {
            rate_c12c12(state.tf, 1.0_rt, rates.fr, rates.rr);
        }
        else if constexpr (rate == C12_O16_to_Mg24_He4 || rate == C12_O16_to_Si28) {
            rate_c12o16(state.tf, 1.0_rt, rates.fr, rates.rr);
        }
        else if constexpr (rate == O16_O16_to_Si28_He4) {
            rate_o16o16(state.tf, 1.0_rt, rates.fr, rates.rr);
        }
        else if constexpr (rate == O16_He4_to_Ne20) {
            rate_o16ag(state.tf, 1.0_rt, rates.fr, rates.rr);
        }
        else if constexpr (rate == Ne20_He4_to_Mg24) {
            rate_ne20ag(state.tf, 1.0_rt, rates.fr, rates.rr);
        }
        else if constexpr (rate == Mg24_He4_to_Si28) {
            rate_mg24ag(state.tf, 1.0_rt, rates.fr, rates.rr);
        }
        else if constexpr (rate == Ca40_He4_to_Ti44) {
            rate_ca40ag(state.tf, 1.0_rt, rates.fr, rates.rr);
        }
    }

    template<int rate, typename number_t>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void postprocess_rate (const rhs_state_t<number_t>& state, rate_t<number_t>& rates,
                           rate_t<number_t>& rates1, [[maybe_unused]] rate_t<number_t>& rates2,
                           [[maybe_unused]] rate_t<number_t>& rates3)
    {
        using namespace Species;
        using namespace Rates;

        if constexpr (rate == Si28_7He4_to_Ni56) {
            // first rate corresponds to Ca40_He4_to_Ti44

            if (static_cast<amrex::Real>(state.tf.t9) > 2.5_rt && state.y(C12) + state.y(O16) <= 4.0e-3_rt) {

                number_t t992  = state.tf.t972 * state.tf.t9;
                number_t t9i92 = 1.0_rt / t992;

                number_t yeff_ca40   = t9i92 * admath::exp(239.42_rt * state.tf.t9i - 74.741_rt);

                number_t yeff_ti44   = t992  * admath::exp(-274.12_rt * state.tf.t9i + 74.914_rt);

                amrex::Real denom = amrex::Math::powi<3>(state.rho * state.y(He4));

                rates.fr   = yeff_ca40 * denom * rates1.fr;

                amrex::Real zz = 1.0_rt / denom;
                rates.rr = admath::min(1.0e10_rt, yeff_ti44 * rates1.rr * zz);
            }
        }
    }

    template<int spec, typename number_t>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    amrex::Real ener_gener_rate ([[maybe_unused]] const rhs_state_t<number_t>& rhs_state, amrex::Real const& dydt)
    {
        return dydt * network::mion<spec>() * C::Legacy::enuc_conv2;
    }

} // namespace RHS

#endif
