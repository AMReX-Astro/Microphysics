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

using namespace network_rp;

AMREX_INLINE
void actual_network_init() {}

constexpr std::string_view network_name = "aprox13";

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
        else if constexpr (spec == S32) {
            return 271.78250e0_rt;
        }
        else if constexpr (spec == Ar36) {
            return 306.72020e0_rt;
        }
        else if constexpr (spec == Ca40) {
            return 342.05680e0_rt;
        }
        else if constexpr (spec == Ti44) {
            return 375.47720e0_rt;
        }
        else if constexpr (spec == Cr48) {
            return 411.46900e0_rt;
        }
        else if constexpr (spec == Fe52) {
            return 447.70800e0_rt;
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
    constexpr int P_index= -1;
    constexpr int H1_index = -1;
    constexpr int N_index = -1;
    constexpr int He4_index = 0;
}
#endif

namespace Rates {
    enum NetworkRates: std::uint8_t {
        He4_He4_He4_to_C12 = 1,
        C12_He4_to_O16,
        C12_C12_to_Ne20_He4,
        C12_O16_to_Mg24_He4,
        C12_O16_to_Si28,
        O16_O16_to_Si28_He4,
        O16_O16_to_P31_P,
        O16_O16_to_S32,
        O16_He4_to_Ne20,
        Ne20_He4_to_Mg24,
        Mg24_He4_to_Si28,
        Si28_He4_to_S32,
        S32_He4_to_Ar36,
        Ar36_He4_to_Ca40,
        Ca40_He4_to_Ti44,
        Ti44_He4_to_Cr48,
        Cr48_He4_to_Fe52,
        Fe52_He4_to_Ni56,
        Mg24_He4_to_Al27_P,
        Al27_P_to_Si28,
        Si28_He4_to_P31_P,
        P31_P_to_S32,
        S32_He4_to_Cl35_P,
        Cl35_P_to_Ar36,
        Ar36_He4_to_K39_P,
        K39_P_to_Ca40,
        Ca40_He4_to_Sc43_P,
        Sc43_P_to_Ti44,
        Ti44_He4_to_V47_P,
        V47_P_to_Cr48,
        Cr48_He4_to_Mn51_P,
        Mn51_P_to_Fe52,
        Fe52_He4_to_Co55_P,
        Co55_P_to_Ni56,
        NumRates=Co55_P_to_Ni56
    };
}

namespace RHS {

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    constexpr rhs_t rhs_data (int rate)
    {
        using namespace Species;
        using namespace Rates;

        rhs_t data;

        switch (rate) {  // NOLINT(bugprone-switch-missing-default-case)

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

            data.forward_branching_ratio = 0.5_rt;
            break;

        case O16_O16_to_Si28_He4:

            // this models the flows:
            //  O16(O16,alpha)Si28
            //  O16(O16,p)P31(p,alpha)Si28

            data.species_A = O16;
            data.species_D = Si28;
            data.species_E = He4;

            data.number_A = 2.0_rt;
            data.number_D = 1.0_rt;
            data.number_E = 1.0_rt;

            data.exponent_A = 2;
            data.exponent_D = 1;
            data.exponent_E = 1;

            data.forward_branching_ratio = 0.56_rt;

            data.additional_reaction_1 = Si28_He4_to_P31_P;
            data.additional_reaction_2 = P31_P_to_S32;
            data.additional_reaction_3 = O16_O16_to_P31_P;
            break;

        case O16_O16_to_P31_P:
            data.species_A = O16;
            data.species_D = P31;
            data.species_E = P;

            data.number_A = 2.0_rt;
            data.number_D = 1.0_rt;
            data.number_E = 1.0_rt;

            data.exponent_A = 2;
            data.exponent_D = 1;
            data.exponent_E = 1;

            data.forward_branching_ratio = 0.34_rt;
            break;

        case O16_O16_to_S32:

            // this models
            //  O16(O16,gamma)S32
            //  O16(O16,p)P31(p,g)S32
            //
            // we also need Si28(alpha,p)P31 to find the equilibrium
            // number of protons

            data.species_A = O16;
            data.species_D = S32;

            data.number_A = 2.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 2;
            data.exponent_D = 1;

            data.forward_branching_ratio = 0.10_rt;

            data.additional_reaction_1 = Si28_He4_to_P31_P;
            data.additional_reaction_2 = P31_P_to_S32;
            data.additional_reaction_3 = O16_O16_to_P31_P;
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

            // this models
            //  Mg24(alpha,gamma)Si28
            //  Mg24(alpha,p)Al27(p,gamma)Si28

            data.species_A = Mg24;
            data.species_B = He4;
            data.species_D = Si28;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;

            data.additional_reaction_1 = Mg24_He4_to_Al27_P;
            data.additional_reaction_2 = Al27_P_to_Si28;
            break;

        case Si28_He4_to_S32:

            // this models
            //  Si28(alpha,gamma)S32
            //  Si28(alpha,p)P31(p,gamma)S32

            data.species_A = Si28;
            data.species_B = He4;
            data.species_D = S32;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;

            data.additional_reaction_1 = Si28_He4_to_P31_P;
            data.additional_reaction_2 = P31_P_to_S32;
            break;

        case S32_He4_to_Ar36:

            // this models
            //  S32(alpha,gamma)Ar36
            //  S32(alpha,p)Cl35(p,gamma)Ar36

            data.species_A = S32;
            data.species_B = He4;
            data.species_D = Ar36;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;

            data.additional_reaction_1 = S32_He4_to_Cl35_P;
            data.additional_reaction_2 = Cl35_P_to_Ar36;
            break;

        case Ar36_He4_to_Ca40:

            // this models
            //  Ar36(alpha,gamma)Ca40
            //  Ar36(alpha,p)K39(p,gamma)Ca40

            data.species_A = Ar36;
            data.species_B = He4;
            data.species_D = Ca40;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;

            data.additional_reaction_1 = Ar36_He4_to_K39_P;
            data.additional_reaction_2 = K39_P_to_Ca40;
            break;

        case Ca40_He4_to_Ti44:

            // this models
            //  Ca40(alpha,gamma)Ti44
            //  Ca40(alpha,p)Sc43(p,gamma)Ti44

            data.species_A = Ca40;
            data.species_B = He4;
            data.species_D = Ti44;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;

            data.additional_reaction_1 = Ca40_He4_to_Sc43_P;
            data.additional_reaction_2 = Sc43_P_to_Ti44;
            break;

        case Ti44_He4_to_Cr48:

            // this models
            //  Ti44(alpha,gamma)Cr48
            //  Ti44(alpha,p)V47(p,gamma)Cr48

            data.species_A = Ti44;
            data.species_B = He4;
            data.species_D = Cr48;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;

            data.additional_reaction_1 = Ti44_He4_to_V47_P;
            data.additional_reaction_2 = V47_P_to_Cr48;
            break;

        case Cr48_He4_to_Fe52:

            // this models
            //  Cr48(alpha,gamma)Fe52
            //  Cr48(alpha,p)Mn51(p,gamma)Fe52

            data.species_A = Cr48;
            data.species_B = He4;
            data.species_D = Fe52;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;

            data.additional_reaction_1 = Cr48_He4_to_Mn51_P;
            data.additional_reaction_2 = Mn51_P_to_Fe52;
            break;

        case Fe52_He4_to_Ni56:

            // this models
            //  Fe52(alpha,gamma)Ni56
            //  Fe52(alpha,p)Co55(p,gamma)Ni56

            data.species_A = Fe52;
            data.species_B = He4;
            data.species_D = Ni56;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;

            data.additional_reaction_1 = Fe52_He4_to_Co55_P;
            data.additional_reaction_2 = Co55_P_to_Ni56;
            break;

        case Mg24_He4_to_Al27_P:
            data.species_A = Mg24;
            data.species_B = He4;
            data.species_D = Al27;
            data.species_E = P;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;
            data.number_E = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            data.exponent_E = 1;
            break;

        case Al27_P_to_Si28:
            data.species_A = Al27;
            data.species_B = P;
            data.species_D = Si28;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            break;

        case Si28_He4_to_P31_P:
            data.species_A = Si28;
            data.species_B = He4;
            data.species_D = P31;
            data.species_E = P;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;
            data.number_E = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            data.exponent_E = 1;
            break;

        case P31_P_to_S32:
            data.species_A = P31;
            data.species_B = P;
            data.species_D = S32;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            break;

        case S32_He4_to_Cl35_P:
            data.species_A = S32;
            data.species_B = He4;
            data.species_D = Cl35;
            data.species_E = P;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;
            data.number_E = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            data.exponent_E = 1;
            break;

        case Cl35_P_to_Ar36:
            data.species_A = Cl35;
            data.species_B = P;
            data.species_D = Ar36;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            break;

        case Ar36_He4_to_K39_P:
            data.species_A = Ar36;
            data.species_B = He4;
            data.species_D = K39;
            data.species_E = P;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;
            data.number_E = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            data.exponent_E = 1;
            break;

        case K39_P_to_Ca40:
            data.species_A = K39;
            data.species_B = P;
            data.species_D = Ca40;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            break;

        case Ca40_He4_to_Sc43_P:
            data.species_A = Ca40;
            data.species_B = He4;
            data.species_D = Sc43;
            data.species_E = P;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;
            data.number_E = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            data.exponent_E = 1;
            break;

        case Sc43_P_to_Ti44:
            data.species_A = Sc43;
            data.species_B = P;
            data.species_D = Ti44;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            break;

        case Ti44_He4_to_V47_P:
            data.species_A = Ti44;
            data.species_B = He4;
            data.species_D = V47;
            data.species_E = P;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;
            data.number_E = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            data.exponent_E = 1;
            break;

        case V47_P_to_Cr48:
            data.species_A = V47;
            data.species_B = P;
            data.species_D = Cr48;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            break;

        case Cr48_He4_to_Mn51_P:
            data.species_A = Cr48;
            data.species_B = He4;
            data.species_D = Mn51;
            data.species_E = P;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;
            data.number_E = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            data.exponent_E = 1;
            break;

        case Mn51_P_to_Fe52:
            data.species_A = Mn51;
            data.species_B = P;
            data.species_D = Fe52;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            break;

        case Fe52_He4_to_Co55_P:
            data.species_A = Fe52;
            data.species_B = He4;
            data.species_D = Co55;
            data.species_E = P;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;
            data.number_E = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            data.exponent_E = 1;
            break;

        case Co55_P_to_Ni56:
            data.species_A = Co55;
            data.species_B = P;
            data.species_D = Ni56;

            data.number_A = 1.0_rt;
            data.number_B = 1.0_rt;
            data.number_D = 1.0_rt;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
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
        else if constexpr (rate == O16_O16_to_Si28_He4 || rate == O16_O16_to_P31_P || rate == O16_O16_to_S32) {
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
        else if constexpr (rate == Si28_He4_to_S32) {
            rate_si28ag(state.tf, 1.0_rt, rates.fr, rates.rr);
        }
        else if constexpr (rate == S32_He4_to_Ar36) {
            rate_s32ag(state.tf, 1.0_rt, rates.fr, rates.rr);
        }
        else if constexpr (rate == Ar36_He4_to_Ca40) {
            rate_ar36ag(state.tf, 1.0_rt, rates.fr, rates.rr);
        }
        else if constexpr (rate == Ca40_He4_to_Ti44) {
            rate_ca40ag(state.tf, 1.0_rt, rates.fr, rates.rr);
        }
        else if constexpr (rate == Ti44_He4_to_Cr48) {
            rate_ti44ag(state.tf, 1.0_rt, rates.fr, rates.rr);
        }
        else if constexpr (rate == Cr48_He4_to_Fe52) {
            rate_cr48ag(state.tf, 1.0_rt, rates.fr, rates.rr);
        }
        else if constexpr (rate == Fe52_He4_to_Ni56) {
            rate_fe52ag(state.tf, 1.0_rt, rates.fr, rates.rr);
        }
        else if constexpr (rate == Mg24_He4_to_Al27_P) {
            rate_mg24ap(state.tf, 1.0_rt, rates.fr, rates.rr);
        }
        else if constexpr (rate == Al27_P_to_Si28) {
            rate_al27pg(state.tf, 1.0_rt, rates.fr, rates.rr);
        }
        else if constexpr (rate == Si28_He4_to_P31_P) {
            rate_si28ap(state.tf, 1.0_rt, rates.fr, rates.rr);
        }
        else if constexpr (rate == P31_P_to_S32) {
            rate_p31pg(state.tf, 1.0_rt, rates.fr, rates.rr);
        }
        else if constexpr (rate == S32_He4_to_Cl35_P) {
            rate_s32ap(state.tf, 1.0_rt, rates.fr, rates.rr);
        }
        else if constexpr (rate == Cl35_P_to_Ar36) {
            rate_cl35pg(state.tf, 1.0_rt, rates.fr, rates.rr);
        }
        else if constexpr (rate == Ar36_He4_to_K39_P) {
            rate_ar36ap(state.tf, 1.0_rt, rates.fr, rates.rr);
        }
        else if constexpr (rate == K39_P_to_Ca40) {
            rate_k39pg(state.tf, 1.0_rt, rates.fr, rates.rr);
        }
        else if constexpr (rate == Ca40_He4_to_Sc43_P) {
            rate_ca40ap(state.tf, 1.0_rt, rates.fr, rates.rr);
        }
        else if constexpr (rate == Sc43_P_to_Ti44) {
            rate_sc43pg(state.tf, 1.0_rt, rates.fr, rates.rr);
        }
        else if constexpr (rate == Ti44_He4_to_V47_P) {
            rate_ti44ap(state.tf, 1.0_rt, rates.fr, rates.rr);
        }
        else if constexpr (rate == V47_P_to_Cr48) {
            rate_v47pg(state.tf, 1.0_rt, rates.fr, rates.rr);
        }
        else if constexpr (rate == Cr48_He4_to_Mn51_P) {
            rate_cr48ap(state.tf, 1.0_rt, rates.fr, rates.rr);
        }
        else if constexpr (rate == Mn51_P_to_Fe52) {
            rate_mn51pg(state.tf, 1.0_rt, rates.fr, rates.rr);
        }
        else if constexpr (rate == Fe52_He4_to_Co55_P) {
            rate_fe52ap(state.tf, 1.0_rt, rates.fr, rates.rr);
        }
        else if constexpr (rate == Co55_P_to_Ni56) {
            rate_co55pg(state.tf, 1.0_rt, rates.fr, rates.rr);
        }
    }

    template<int rate, typename number_t>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void postprocess_rate ([[maybe_unused]] const rhs_state_t<number_t>& state, rate_t<number_t>& rates,
                           rate_t<number_t>& rates1, rate_t<number_t>& rates2,
                           rate_t<number_t>& rates3)
    {
        using namespace Species;
        using namespace Rates;

        // Each of the following reactions has both an (a,g) path and an
        // equilibrium sequence using (a,p),(p,g). The first path has already
        // been handled and we need to add the second now. For example, for
        // reaction Mg24_He4_to_Si28, the first additional rate is Mg24_He4_to_Al27_P
        // and the second is Al27_P_to_Si28. The algorithm is that we compute
        // the effective equilibrium ratio of the forward and reverse rates for
        // the equilibrium sequence and then add the corresponding fraction of
        // those equilibrium chains to the forward and reverse terms for the
        // main alpha chain term.

        if constexpr (rate == Mg24_He4_to_Si28 ||
                      rate == Si28_He4_to_S32 ||
                      rate == S32_He4_to_Ar36 ||
                      rate == Ar36_He4_to_Ca40 ||
                      rate == Ca40_He4_to_Ti44 ||
                      rate == Ti44_He4_to_Cr48 ||
                      rate == Cr48_He4_to_Fe52 ||
                      rate == Fe52_He4_to_Ni56)
        {
            number_t rate_ir    = 0.0e0_rt;
            number_t denom      = rates1.rr + rates2.fr;

            if (denom > 1.0e-30_rt) {
                number_t zz = 1.0_rt / denom;
                rate_ir    = rates1.rr * zz;
            }

            rates.fr += rates1.fr * (1.0_rt - rate_ir);
            rates.rr += rates2.rr * rate_ir;
        }

        // The O16+O16->Si28+He4 reaction has an additional contribution from the O16+O16->P31+P
        // reaction followed by the reverse Si28+He4<->P31+P reaction.
        if constexpr (rate == O16_O16_to_Si28_He4)
        {
            number_t rate_ir    = 0.0e0_rt;
            number_t denom      = rates1.rr + rates2.fr;

            if (denom > 1.0e-30_rt) {
                number_t zz = 1.0_rt / denom;
                rate_ir    = rates1.rr * zz;
            }

            rates.fr += rates3.fr * rate_ir;
        }

        // The O16+O16->S32 reaction has an additional contribution from the O16+O16->P31+P
        // reaction followed by the forward P31+P<->S32 reaction.
        if constexpr (rate == O16_O16_to_S32)
        {
            number_t rate_ir    = 0.0e0_rt;
            number_t denom      = rates1.rr + rates2.fr;

            if (denom > 1.0e-30_rt) {
                number_t zz = 1.0_rt / denom;
                rate_ir    = rates1.rr * zz;
            }

            rates.fr += rates3.fr * (1.0_rt - rate_ir);
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
