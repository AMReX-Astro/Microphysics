#ifndef actual_network_H
#define actual_network_H

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

const std::string network_name = "aprox21";

namespace network
{
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> bion;
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> mion;
}


namespace Rates {
    enum NetworkRates {
        H1_H1_to_He3 = 1,
        H1_H1_H1_to_He3,
        P_to_N,
        P_N_to_H2,
        P_P_N_N_to_He4,
        H2_P_to_He3,
        He3_N_to_He4,
        He3_He3_to_He4_2H1,
        He3_He4_H1_to_2He4,
        He4_He4_He4_to_C12,
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
        Fe52_N_to_Fe53,
        Fe53_N_to_Fe54,
        Ni56_to_Fe56,
        Fe56_to_Cr56,
        C12_H1_to_N13,
        C12_2H1_to_N14,
        N14_H1_to_O15,
        N15_H1_to_O16,
        N15_H1_to_C12_He4,
        N14_2H1_to_O16,
        N14_2H1_to_C12_He4,
        O16_H1_to_F17,
        O16_2H1_to_N14_He4,
        N14_He4_to_F18,
        N14_He4_to_Ne20,
        Fe52_2N_to_Fe54,
        Fe52_He4_to_Fe54_2P,
        Fe54_2P_to_Ni56,
        Fe54_P_to_Co55,
        Fe54_N_to_Fe55,
        Fe55_N_to_Fe56,
        Fe54_2N_to_Fe56,
        Fe54_He4_to_Co57_P,
        Fe56_P_to_Co57,
        Fe54_He4_to_Fe56_2P,
        NumRates=Fe54_He4_to_Fe56_2P
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

        case P_to_N:
            // irpen and irnep in the original aprox21
            data.species_A = P;
            data.species_D = N;

            data.number_A = 1;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_D = 1;

            data.screen_forward_reaction = 0;
            data.screen_reverse_reaction = 0;

            data.rate_can_be_tabulated = 0;
            break;

        case P_N_to_H2:
            // irhng and irdgn in the original aprox21
            // Represents p(n,g)H2
            // In this network this is entirely an intermediate
            // reaction and does not directly enter into the RHS
            // for P or N. We handle this by setting the
            // number of species consumed and produced to zero.
            data.species_A = P;
            data.species_B = N;
            data.species_D = H2;

            data.number_A = 0;
            data.number_B = 0;
            data.number_D = 0;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;

            data.screen_forward_reaction = 0;
            data.screen_reverse_reaction = 0;
            break;

        case P_P_N_N_to_He4:
            // iralf1 and iralf2 in the original aprox21
            // Represents p(n,g)H2(n,g)H3(p,g)He4 and p(n,g)H2(p,g)He3(n,g)He4
            data.species_A = P;
            data.species_B = N;
            data.species_D = He4;

            data.number_A = 2;
            data.number_B = 2;
            data.number_D = 1;

            data.exponent_A = 2;
            data.exponent_B = 2;
            data.exponent_D = 1;

            data.screen_forward_reaction = 0;
            data.screen_reverse_reaction = 0;

            data.rate_can_be_tabulated = 0;

            data.apply_identical_particle_factor = 0;

            data.additional_reaction_1 = P_N_to_H2;
            data.additional_reaction_2 = H2_P_to_He3;
            data.additional_reaction_3 = He3_N_to_He4;
            break;

        case H1_H1_to_He3:
            // irpp in the original aprox21
            // Represents the standard proton-proton chain,
            // H1(H1,e+nu)H2(H1,g)He3
            data.species_A = H1;
            data.species_D = He3;

            data.number_A = 3;
            data.number_D = 1;

            data.exponent_A = 2;
            data.exponent_D = 1;
            break;

        case H1_H1_H1_to_He3:
            // Represents H1(e+H1)H2(H1,g)He3
            // This is a variant of the above reaction; the end
            // result is still three protons to He3, but because
            // the proton is first converted to a neutron, the
            // identical particle factor would not apply. In the
            // original aprox21, this is one of the two reactions
            // that uses the electron capture rate irpen.
            data.species_A = H1;
            data.species_D = He3;

            data.number_A = 3;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_D = 1;

            data.rate_can_be_tabulated = 0;

            data.screen_forward_reaction = 0;
            data.screen_reverse_reaction = 0;

            data.additional_reaction_1 = P_to_N;
            break;

        case H2_P_to_He3:
            // irdpg and irhegp in the original aprox21
            // Represents H2(p,g)He3
            // In this network this is entirely an intermediate
            // reaction and does not directly enter into the RHS
            // for P or He3. We handle this by setting the
            // number of species consumed and produced to zero.
            data.species_A = H2;
            data.species_B = P;
            data.species_D = He3;

            data.number_A = 0;
            data.number_B = 0;
            data.number_D = 0;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            break;

        case He3_N_to_He4:
            // irheng and irhegn in the original aprox21
            // Represents He3(n,g)He4
            // In this network this is entirely an intermediate
            // reaction and does not directly enter into the RHS
            // for He3, N, or He4. We handle this by setting the
            // number of species consumed and produced to zero.
            data.species_A = He3;
            data.species_B = N;
            data.species_D = He4;

            data.number_A = 0;
            data.number_B = 0;
            data.number_D = 0;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;

            data.screen_forward_reaction = 0;
            data.screen_reverse_reaction = 0;
            break;

        case He3_He3_to_He4_2H1:
            // ir33 in the original aprox21
            data.species_A = He3;
            data.species_D = He4;
            data.species_E = H1;

            data.number_A = 2;
            data.number_D = 1;
            data.number_E = 2;

            data.exponent_A = 2;
            data.exponent_D = 1;
            data.exponent_E = 1;
            break;

        case He3_He4_H1_to_2He4:
            // irhe3ag in the original aprox21
            // Represents He3(a,g)Be7(p,g)B8(e+nu)Be8(2a)
            data.species_A = He3;
            data.species_B = He4;
            data.species_C = H1;
            data.species_D = He4;

            data.number_A = 1;
            data.number_B = 1;
            data.number_C = 1;
            data.number_D = 2;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_C = 0;
            data.exponent_D = 1;

            data.apply_identical_particle_factor = 0;
            break;

        case He4_He4_He4_to_C12:
            data.species_A = He4;
            data.species_D = C12;

            data.number_A = 3;
            data.number_D = 1;

            data.exponent_A = 3;
            data.exponent_D = 1;
            break;

        case C12_H1_to_N13:
            // ircpg in the original aprox21
            data.species_A = C12;
            data.species_B = H1;
            data.species_D = N13;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            break;

        case C12_2H1_to_N14:
            // Represents C12(p,g)N13(e+nu)C13(p,g)N14
            data.species_A = C12;
            data.species_B = H1;
            data.species_D = N14;

            data.number_A = 1;
            data.number_B = 2;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;

            data.screen_forward_reaction = 0;
            data.screen_reverse_reaction = 0;

            data.rate_can_be_tabulated = 0;

            data.additional_reaction_1 = C12_H1_to_N13;
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

            data.screen_forward_reaction = 0;
            data.screen_reverse_reaction = 0;

            data.rate_can_be_tabulated = 0;

            data.additional_reaction_1 = C12_O16_to_Mg24_He4;
            break;

        case N14_H1_to_O15:
            // irnpg in the original aprox21
            data.species_A = N14;
            data.species_B = H1;
            data.species_D = O15;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            break;

        case N15_H1_to_O16:
            // ifg in the original aprox21
            data.species_A = N15;
            data.species_B = H1;
            data.species_D = O16;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;

            data.screen_forward_reaction = 0;
            data.screen_reverse_reaction = 0;
            break;

        case N15_H1_to_C12_He4:
            // ifa in the original aprox21
            data.species_A = N15;
            data.species_B = H1;
            data.species_D = C12;
            data.species_E = He4;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;
            data.number_E = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            data.exponent_E = 1;

            data.screen_forward_reaction = 0;
            data.screen_reverse_reaction = 0;
            break;

        case N14_2H1_to_O16:
            // Represents N14(p,g)O15(e+nu)N15(p,g)O16
            data.species_A = N14;
            data.species_B = H1;
            data.species_D = O16;

            data.number_A = 1;
            data.number_B = 2;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;

            data.screen_forward_reaction = 0;
            data.screen_reverse_reaction = 0;

            data.rate_can_be_tabulated = 0;

            data.additional_reaction_1 = N14_H1_to_O15;
            data.additional_reaction_2 = N15_H1_to_O16;
            data.additional_reaction_3 = N15_H1_to_C12_He4;
            break;

        case N14_2H1_to_C12_He4:
            // Represents N14(p,g)O15(e+nu)N15(p,a)C12
            data.species_A = N14;
            data.species_B = H1;
            data.species_D = C12;
            data.species_E = He4;

            data.number_A = 1;
            data.number_B = 2;
            data.number_D = 1;
            data.number_E = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            data.exponent_E = 1;

            data.screen_forward_reaction = 0;
            data.screen_reverse_reaction = 0;

            data.rate_can_be_tabulated = 0;

            data.additional_reaction_1 = N14_H1_to_O15;
            data.additional_reaction_2 = N15_H1_to_O16;
            data.additional_reaction_3 = N15_H1_to_C12_He4;
            break;

        case N14_He4_to_F18:
            // irnag in the original aprox21
            // Represents N14 + He4 -> F18
            data.species_A = N14;
            data.species_B = He4;
            data.species_D = F18;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            break;

        case N14_He4_to_Ne20:
            // Represents N14 + 1.5 He4 -> Ne20
            // or, equivalently, 2 N14 + 3 He4 -> 2 Ne20
            data.species_A = N14;
            data.species_B = He4;
            data.species_D = Ne20;

            data.number_A = 2;
            data.number_B = 3;
            data.number_D = 2;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;

            data.screen_forward_reaction = 0;
            data.screen_reverse_reaction = 0;

            data.rate_can_be_tabulated = 0;

            data.additional_reaction_1 = N14_He4_to_F18;
            break;

        case O16_H1_to_F17:
            // iropg in the original aprox21
            data.species_A = O16;
            data.species_B = H1;
            data.species_D = F17;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            break;

        case O16_2H1_to_N14_He4:
            // Represents O16(p,g)F17(e+nu)O17(p,a)N14
            data.species_A = O16;
            data.species_B = H1;
            data.species_D = N14;
            data.species_E = He4;

            data.number_A = 1;
            data.number_B = 2;
            data.number_D = 1;
            data.number_E = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            data.exponent_E = 1;

            data.screen_forward_reaction = 0;
            data.screen_reverse_reaction = 0;

            data.rate_can_be_tabulated = 0;

            data.additional_reaction_1 = O16_H1_to_F17;
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

            data.screen_forward_reaction = 0;
            data.screen_reverse_reaction = 0;

            data.rate_can_be_tabulated = 0;

            data.additional_reaction_1 = Si28_He4_to_P31_P;
            data.additional_reaction_2 = P31_P_to_S32;
            data.additional_reaction_3 = O16_O16_to_P31_P;
            break;

        case O16_O16_to_P31_P:
            data.species_A = O16;
            data.species_D = P31;
            data.species_E = P;

            data.number_A = 2;
            data.number_D = 1;
            data.number_E = 1;

            data.exponent_A = 2;
            data.exponent_D = 1;
            data.exponent_E = 1;
            break;

        case O16_O16_to_S32:
            data.species_A = O16;
            data.species_D = S32;

            data.number_A = 2;
            data.number_D = 1;

            data.exponent_A = 2;
            data.exponent_D = 1;

            data.screen_forward_reaction = 0;
            data.screen_reverse_reaction = 0;

            data.rate_can_be_tabulated = 0;

            data.additional_reaction_1 = Si28_He4_to_P31_P;
            data.additional_reaction_2 = P31_P_to_S32;
            data.additional_reaction_3 = O16_O16_to_P31_P;
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

            data.additional_reaction_1 = Mg24_He4_to_Al27_P;
            data.additional_reaction_2 = Al27_P_to_Si28;
            break;

        case Si28_He4_to_S32:
            data.species_A = Si28;
            data.species_B = He4;
            data.species_D = S32;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;

            data.additional_reaction_1 = Si28_He4_to_P31_P;
            data.additional_reaction_2 = P31_P_to_S32;
            break;

        case S32_He4_to_Ar36:
            data.species_A = S32;
            data.species_B = He4;
            data.species_D = Ar36;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;

            data.additional_reaction_1 = S32_He4_to_Cl35_P;
            data.additional_reaction_2 = Cl35_P_to_Ar36;
            break;

        case Ar36_He4_to_Ca40:
            data.species_A = Ar36;
            data.species_B = He4;
            data.species_D = Ca40;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;

            data.additional_reaction_1 = Ar36_He4_to_K39_P;
            data.additional_reaction_2 = K39_P_to_Ca40;
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

            data.additional_reaction_1 = Ca40_He4_to_Sc43_P;
            data.additional_reaction_2 = Sc43_P_to_Ti44;
            break;

        case Ti44_He4_to_Cr48:
            data.species_A = Ti44;
            data.species_B = He4;
            data.species_D = Cr48;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;

            data.additional_reaction_1 = Ti44_He4_to_V47_P;
            data.additional_reaction_2 = V47_P_to_Cr48;
            break;

        case Cr48_He4_to_Fe52:
            data.species_A = Cr48;
            data.species_B = He4;
            data.species_D = Fe52;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;

            data.additional_reaction_1 = Cr48_He4_to_Mn51_P;
            data.additional_reaction_2 = Mn51_P_to_Fe52;
            break;

        case Fe52_He4_to_Ni56:
            data.species_A = Fe52;
            data.species_B = He4;
            data.species_D = Ni56;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;

            data.additional_reaction_1 = Fe52_He4_to_Co55_P;
            data.additional_reaction_2 = Co55_P_to_Ni56;
            data.additional_reaction_3 = Fe54_P_to_Co55;
            break;

        case Mg24_He4_to_Al27_P:
            data.species_A = Mg24;
            data.species_B = He4;
            data.species_D = Al27;
            data.species_E = P;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;
            data.number_E = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            data.exponent_E = 1;
            break;

        case Al27_P_to_Si28:
            data.species_A = Al27;
            data.species_B = P;
            data.species_D = Si28;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            break;

        case Si28_He4_to_P31_P:
            data.species_A = Si28;
            data.species_B = He4;
            data.species_D = P31;
            data.species_E = P;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;
            data.number_E = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            data.exponent_E = 1;
            break;

        case P31_P_to_S32:
            data.species_A = P31;
            data.species_B = P;
            data.species_D = S32;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            break;

        case S32_He4_to_Cl35_P:
            data.species_A = S32;
            data.species_B = He4;
            data.species_D = Cl35;
            data.species_E = P;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;
            data.number_E = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            data.exponent_E = 1;
            break;

        case Cl35_P_to_Ar36:
            data.species_A = Cl35;
            data.species_B = P;
            data.species_D = Ar36;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            break;

        case Ar36_He4_to_K39_P:
            data.species_A = Ar36;
            data.species_B = He4;
            data.species_D = K39;
            data.species_E = P;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;
            data.number_E = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            data.exponent_E = 1;
            break;

        case K39_P_to_Ca40:
            data.species_A = K39;
            data.species_B = P;
            data.species_D = Ca40;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            break;

        case Ca40_He4_to_Sc43_P:
            data.species_A = Ca40;
            data.species_B = He4;
            data.species_D = Sc43;
            data.species_E = P;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;
            data.number_E = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            data.exponent_E = 1;
            break;

        case Sc43_P_to_Ti44:
            data.species_A = Sc43;
            data.species_B = P;
            data.species_D = Ti44;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            break;

        case Ti44_He4_to_V47_P:
            data.species_A = Ti44;
            data.species_B = He4;
            data.species_D = V47;
            data.species_E = P;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;
            data.number_E = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            data.exponent_E = 1;
            break;

        case V47_P_to_Cr48:
            data.species_A = V47;
            data.species_B = P;
            data.species_D = Cr48;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            break;

        case Cr48_He4_to_Mn51_P:
            data.species_A = Cr48;
            data.species_B = He4;
            data.species_D = Mn51;
            data.species_E = P;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;
            data.number_E = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            data.exponent_E = 1;
            break;

        case Mn51_P_to_Fe52:
            data.species_A = Mn51;
            data.species_B = P;
            data.species_D = Fe52;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            break;

        case Fe52_He4_to_Co55_P:
            data.species_A = Fe52;
            data.species_B = He4;
            data.species_D = Co55;
            data.species_E = P;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;
            data.number_E = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            data.exponent_E = 1;
            break;

        case Co55_P_to_Ni56:
            data.species_A = Co55;
            data.species_B = P;
            data.species_D = Ni56;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            break;

        case Fe52_N_to_Fe53:
            data.species_A = Fe52;
            data.species_B = N;
            data.species_D = Fe53;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;

            data.screen_forward_reaction = 0;
            data.screen_reverse_reaction = 0;
            break;

        case Fe53_N_to_Fe54:
            data.species_A = Fe53;
            data.species_B = N;
            data.species_D = Fe54;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;

            data.screen_forward_reaction = 0;
            data.screen_reverse_reaction = 0;
            break;

        case Fe52_2N_to_Fe54:
            data.species_A = Fe52;
            data.species_B = N;
            data.species_D = Fe54;

            data.number_A = 1;
            data.number_B = 2;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 2;
            data.exponent_D = 1;

            data.screen_forward_reaction = 0;
            data.screen_reverse_reaction = 0;

            data.rate_can_be_tabulated = 0;

            data.apply_identical_particle_factor = 0;

            data.additional_reaction_1 = Fe52_N_to_Fe53;
            data.additional_reaction_2 = Fe53_N_to_Fe54;
            break;

        case Fe54_P_to_Co55:
            data.species_A = Fe54;
            data.species_B = P;
            data.species_D = Co55;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            break;

        case Fe52_He4_to_Fe54_2P:
            data.species_A = Fe52;
            data.species_B = He4;
            data.species_D = Fe54;
            data.species_E = P;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;
            data.number_E = 2;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            data.exponent_E = 2;

            data.screen_forward_reaction = 0;
            data.screen_reverse_reaction = 0;

            data.rate_can_be_tabulated = 0;

            data.apply_identical_particle_factor = 0;

            data.additional_reaction_1 = Fe52_He4_to_Co55_P;
            data.additional_reaction_2 = Fe54_P_to_Co55;
            data.additional_reaction_3 = Co55_P_to_Ni56;
            break;

        case Fe54_2P_to_Ni56:
            data.species_A = Fe54;
            data.species_B = P;
            data.species_D = Ni56;

            data.number_A = 1;
            data.number_B = 2;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 2;
            data.exponent_D = 1;

            data.screen_forward_reaction = 0;
            data.screen_reverse_reaction = 0;

            data.rate_can_be_tabulated = 0;

            data.apply_identical_particle_factor = 0;

            data.additional_reaction_1 = Fe54_P_to_Co55;
            data.additional_reaction_2 = Co55_P_to_Ni56;
            data.additional_reaction_3 = Fe52_He4_to_Co55_P;
            break;

        case Ni56_to_Fe56:
            data.species_A = Ni56;
            data.species_D = Fe56;

            data.number_A = 1;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_D = 1;

            data.rate_can_be_tabulated = 0;

            data.screen_forward_reaction = 0;
            data.screen_reverse_reaction = 0;
            break;

        case Fe56_to_Cr56:
            data.species_A = Fe56;
            data.species_D = Cr56;

            data.number_A = 1;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_D = 1;

            data.rate_can_be_tabulated = 0;

            data.screen_forward_reaction = 0;
            data.screen_reverse_reaction = 0;

            data.additional_reaction_1 = Ni56_to_Fe56;
            break;

        case Fe54_N_to_Fe55:
            data.species_A = Fe54;
            data.species_B = N;
            data.species_D = Fe55;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;

            data.screen_forward_reaction = 0;
            data.screen_reverse_reaction = 0;
            break;

        case Fe55_N_to_Fe56:
            data.species_A = Fe55;
            data.species_B = N;
            data.species_D = Fe56;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;

            data.screen_forward_reaction = 0;
            data.screen_reverse_reaction = 0;
            break;

        case Fe54_2N_to_Fe56:
            data.species_A = Fe54;
            data.species_B = N;
            data.species_D = Fe56;

            data.number_A = 1;
            data.number_B = 2;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 2;
            data.exponent_D = 1;

            data.screen_forward_reaction = 0;
            data.screen_reverse_reaction = 0;

            data.rate_can_be_tabulated = 0;

            data.apply_identical_particle_factor = 0;

            data.additional_reaction_1 = Fe54_N_to_Fe55;
            data.additional_reaction_2 = Fe55_N_to_Fe56;
            break;

        case Fe54_He4_to_Co57_P:
            data.species_A = Fe54;
            data.species_B = He4;
            data.species_D = Co57;
            data.species_E = P;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;
            data.number_E = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            data.exponent_E = 1;
            break;

        case Fe56_P_to_Co57:
            data.species_A = Fe56;
            data.species_B = P;
            data.species_D = Co57;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            break;

        case Fe54_He4_to_Fe56_2P:
            data.species_A = Fe54;
            data.species_B = He4;
            data.species_D = Fe56;
            data.species_E = P;

            data.number_A = 1;
            data.number_B = 1;
            data.number_D = 1;
            data.number_E = 2;

            data.exponent_A = 1;
            data.exponent_B = 1;
            data.exponent_D = 1;
            data.exponent_E = 2;

            data.screen_forward_reaction = 0;
            data.screen_reverse_reaction = 0;

            data.rate_can_be_tabulated = 0;

            data.apply_identical_particle_factor = 0;

            data.additional_reaction_1 = Fe54_He4_to_Co57_P;
            data.additional_reaction_2 = Fe56_P_to_Co57;
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

        if constexpr (rate == H1_H1_to_He3) {
            Real tmp1, tmp2;
            rate_pp(state.tf, 1.0_rt, rates.fr, rates.frdt, tmp1, tmp2);
        }
        else if constexpr (rate == P_to_N) {
            if (state.tf.temp >= 1.0e6_rt && state.rho >= 1.0e-9_rt) {
                Real tmp1, tmp2;
                ecapnuc(state.eta, state.tf.temp, rates.fr, rates.rr, tmp1, tmp2);
            }
        }
        else if constexpr (rate == P_N_to_H2) {
            rate_png(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == H2_P_to_He3) {
            rate_dpg(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == He3_N_to_He4) {
            rate_he3ng(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == He3_He3_to_He4_2H1) {
            Real tmp1, tmp2;
            rate_he3he3(state.tf, 1.0_rt, rates.fr, rates.frdt, tmp1, tmp2);
        }
        else if constexpr (rate == He3_He4_H1_to_2He4) {
            Real tmp1, tmp2;
            rate_he3he4(state.tf, 1.0_rt, rates.fr, rates.frdt, tmp1, tmp2);
        }
        else if constexpr (rate == C12_He4_to_O16) {
            rate_c12ag(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == He4_He4_He4_to_C12) {
            rate_triplealf(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == C12_H1_to_N13) {
            rate_c12pg(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == C12_C12_to_Ne20_He4) {
            rate_c12c12(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == C12_O16_to_Mg24_He4) {
            rate_c12o16(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == N14_H1_to_O15) {
            rate_n14pg(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == N15_H1_to_O16) {
            rate_n15pg(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == N15_H1_to_C12_He4) {
            rate_n15pa(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == N14_He4_to_F18) {
            rate_n14ag(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == O16_H1_to_F17) {
            rate_o16pg(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == O16_O16_to_P31_P) {
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
        else if constexpr (rate == Si28_He4_to_S32) {
            rate_si28ag(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == S32_He4_to_Ar36) {
            rate_s32ag(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == Ar36_He4_to_Ca40) {
            rate_ar36ag(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == Ca40_He4_to_Ti44) {
            rate_ca40ag(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == Ti44_He4_to_Cr48) {
            rate_ti44ag(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == Cr48_He4_to_Fe52) {
            rate_cr48ag(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == Fe52_He4_to_Ni56) {
            rate_fe52ag(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == Mg24_He4_to_Al27_P) {
            rate_mg24ap(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == Al27_P_to_Si28) {
            rate_al27pg(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == Si28_He4_to_P31_P) {
            rate_si28ap(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == P31_P_to_S32) {
            rate_p31pg(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == S32_He4_to_Cl35_P) {
            rate_s32ap(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == Cl35_P_to_Ar36) {
            rate_cl35pg(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == Ar36_He4_to_K39_P) {
            rate_ar36ap(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == K39_P_to_Ca40) {
            rate_k39pg(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == Ca40_He4_to_Sc43_P) {
            rate_ca40ap(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == Sc43_P_to_Ti44) {
            rate_sc43pg(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == Ti44_He4_to_V47_P) {
            rate_ti44ap(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == V47_P_to_Cr48) {
            rate_v47pg(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == Cr48_He4_to_Mn51_P) {
            rate_cr48ap(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == Mn51_P_to_Fe52) {
            rate_mn51pg(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == Fe52_He4_to_Co55_P) {
            rate_fe52ap(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == Co55_P_to_Ni56) {
            rate_co55pg(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == Fe52_N_to_Fe53) {
            rate_fe52ng(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == Fe53_N_to_Fe54) {
            rate_fe53ng(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == Fe54_P_to_Co55) {
            rate_fe54pg(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == Ni56_to_Fe56) {
            if (state.tf.temp >= 1.0e6_rt && state.rho >= 1.0e-9_rt) {
                Real tmp;
                langanke(state.tf.temp, state.rho, state.y(Ni56), state.y_e, rates.fr, tmp);
            }
        }
        else if constexpr (rate == Fe54_N_to_Fe55) {
            rate_fe54ng(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == Fe55_N_to_Fe56) {
            rate_fe55ng(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == Fe54_He4_to_Co57_P) {
            rate_fe54ap(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
        else if constexpr (rate == Fe56_P_to_Co57) {
            rate_fe56pg(state.tf, 1.0_rt, rates.fr, rates.frdt, rates.rr, rates.rrdt);
        }
    }

    template<int rate>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void postprocess_rate (const rhs_state_t& state, rate_t& rates,
                           rate_t& rates1, rate_t& rates2, rate_t& rates3)
    {
        using namespace Species;
        using namespace Rates;

        if constexpr (rate == P_P_N_N_to_He4)
        {
            // In the original aprox21:
            // Rate 1 == P_N_to_H2 == irhng (forward), irdgn (reverse)
            // Rate 2 == H2_P_to_He3 == irdpg (forward), irhegp (reverse)
            // Rate 3 == He3_N_to_He4 == irheng (forward), irhegn (reverse)

            Real denom = rates2.rr * rates1.rr + state.y(N) * rates3.fr * rates1.rr + state.y(N) * state.y(P) * rates3.fr * rates2.fr;

            if (denom > 1.0e-50_rt && state.tf.t9 > 1.5_rt) {
                Real denomdt = rates2.rrdt * rates1.rr + rates2.rr * rates1.rrdt +
                               state.y(N) * (rates3.frdt * rates1.rr + rates3.fr * rates1.rrdt) +
                               state.y(N) * state.y(P) * (rates3.frdt * rates2.fr + rates3.fr * rates2.frdt);

                Real zz = 1.0_rt / denom;

                // iralf1 in the original aprox21
                rates.rr = rates3.rr * rates2.rr * rates1.rr * zz;
                rates.rrdt = rates3.rrdt * rates2.rr * rates1.rr * zz +
                             rates3.rr * rates2.rrdt * rates1.rr * zz +
                             rates3.rr * rates2.rr * rates1.rrdt * zz -
                             rates.rr * zz * denomdt;

                // iralf2 in the original aprox21
                rates.fr = rates3.fr * rates2.fr * rates1.fr * zz;
                rates.frdt = rates3.frdt * rates2.fr * rates1.fr * zz +
                             rates3.fr * rates2.frdt * rates1.fr * zz +
                             rates3.fr * rates2.fr * rates1.frdt * zz -
                             rates.fr * zz * denomdt;
            }
        }

        // Set this equal to the p+e->n rate.

        if constexpr (rate == H1_H1_H1_to_He3) {
            rates.fr = rates1.fr;
        }

        // Beta limit He3 + He4 by the B8 decay half-life

        if constexpr (rate == He3_He4_H1_to_2He4)
        {
            Real xx = 0.896_rt / state.y(He4);
            rates.fr = amrex::min(rates.fr, xx);
            if (rates.fr == xx) {
                rates.frdt = 0.0_rt;
            }
        }

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
                      rate == Cr48_He4_to_Fe52)
        {
            Real rate_ir    = 0.0e0_rt;
            Real dratedt_ir = 0.0e0_rt;
            Real denom      = rates1.rr + rates2.fr;
            Real denomdt    = rates1.rrdt + rates2.frdt;

            if (denom > 1.0e-30_rt) {
                Real zz = 1.0_rt / denom;
                rate_ir    = rates1.rr * zz;
                dratedt_ir = (rates1.rrdt - rate_ir * denomdt) * zz;
            }

            rates.fr += rates1.fr * (1.0_rt - rate_ir);
            rates.frdt += rates1.frdt * (1.0_rt - rate_ir) - rates1.fr * dratedt_ir;

            rates.rr += rates2.rr * rate_ir;
            rates.rrdt += rates2.rrdt * rate_ir + rates2.rr * dratedt_ir;
        }

        // The O16+O16->Si28+He4 reaction has an additional contribution from the O16+O16->P31+P
        // reaction followed by the reverse Si28+He4<->P31+P reaction.
        if constexpr (rate == O16_O16_to_Si28_He4)
        {
            Real rate_ir    = 0.0e0_rt;
            Real dratedt_ir = 0.0e0_rt;
            Real denom      = rates1.rr + rates2.fr;
            Real denomdt    = rates1.rrdt + rates2.frdt;

            if (denom > 1.0e-30_rt) {
                Real zz = 1.0_rt / denom;
                rate_ir    = rates1.rr * zz;
                dratedt_ir = (rates1.rrdt - rate_ir * denomdt) * zz;
            }

            rates.fr = 0.56_rt * rates3.fr + 0.34_rt * rates3.fr * rate_ir;
            rates.frdt = 0.56_rt * rates3.frdt + 0.34_rt * rates3.frdt * rate_ir + 0.34_rt * rates3.fr * dratedt_ir;
        }

        // The O16+O16->S32 reaction has an additional contribution from the O16+O16->P31+P
        // reaction followed by the forward P31+P<->S32 reaction.
        if constexpr (rate == O16_O16_to_S32)
        {
            Real rate_ir    = 0.0e0_rt;
            Real dratedt_ir = 0.0e0_rt;
            Real denom      = rates1.rr + rates2.fr;
            Real denomdt    = rates1.rrdt + rates2.frdt;

            if (denom > 1.0e-30_rt) {
                Real zz = 1.0_rt / denom;
                rate_ir    = rates1.rr * zz;
                dratedt_ir = (rates1.rrdt - rate_ir * denomdt) * zz;
            }

            rates.fr = 0.1_rt * rates3.fr + 0.34_rt * rates3.fr * (1.0_rt - rate_ir);
            rates.frdt = 0.1_rt * rates3.frdt + 0.34_rt * rates3.frdt * (1.0_rt - rate_ir) - 0.34_rt * rates3.fr * dratedt_ir;
        }

        if constexpr (rate == C12_2H1_to_N14)
        {
            rates.fr = rates1.fr;
            rates.frdt = rates1.frdt;
        }

        if constexpr (rate == C12_O16_to_Si28)
        {
            rates.fr = rates1.fr;
            rates.frdt = rates1.frdt;
            rates.rr = rates1.rr;
            rates.rrdt = rates1.rrdt;
        }

        if constexpr (rate == N14_2H1_to_O16)
        {
            // Rate 1: N14_H1_to_O15

            // Beta limit

            if (state.y(H1) > 1.0e-30_rt) {
                Real xx = 5.68e-3_rt / (state.y(H1) * 1.57_rt);
                rates1.fr = amrex::min(rates1.fr, xx);
                if (rates1.fr == xx) {
                    rates1.frdt = 0.0_rt;
                }
            }

            // Rate 2: O15_H1_to_O16
            // Rate 3: O15_H1_to_C12_He4

            Real tot = rates2.fr + rates3.fr;
            Real invtot = 1.0_rt / tot;

            rates.fr = rates1.fr * (rates2.fr * invtot);

            Real dtotdt = rates2.frdt + rates3.frdt;

            rates.frdt = rates1.fr * (rates2.frdt * invtot - rates2.fr * invtot * invtot * dtotdt);
        }

        if constexpr (rate == N14_2H1_to_C12_He4)
        {
            // Rate 1: N14_H1_to_O15

            // Beta limit

            if (state.y(H1) > 1.0e-30_rt) {
                Real xx = 5.68e-3_rt / (state.y(H1) * 1.57_rt);
                rates1.fr = amrex::min(rates1.fr, xx);
                if (rates1.fr == xx) {
                    rates1.frdt = 0.0_rt;
                }
            }

            // Rate 2: O15_H1_to_O16
            // Rate 3: O15_H1_to_C12_He4

            Real tot = rates2.fr + rates3.fr;
            Real invtot = 1.0_rt / tot;

            rates.fr = rates1.fr * (rates3.fr * invtot);

            Real dtotdt = rates2.frdt + rates3.frdt;

            rates.frdt = rates1.fr * (rates3.frdt * invtot - rates3.fr * invtot * invtot * dtotdt);
        }

        if constexpr (rate == N14_He4_to_Ne20)
        {
            // We doubled the number of nucleons per reaction
            // in defining the reaction so that the number of reactants
            // on each side would be integers, so to compensate,
            // we'll halve the rate.

            rates.fr = rates1.fr * 0.5_rt;
            rates.frdt = rates1.frdt * 0.5_rt;
        }

        if constexpr (rate == O16_2H1_to_N14_He4)
        {
            // Beta limit

            Real xx = 0.0105_rt / state.y(H1);
            rates1.fr = amrex::min(rates1.fr, xx);
            if (rates1.fr == xx) {
                rates1.frdt = 0.0_rt;
            }

            rates.fr = rates1.fr;
            rates.frdt = rates1.frdt;
        }

        if constexpr (rate == Fe52_2N_to_Fe54)
        {
            // In the original aprox21:
            // Rate 1 == Fe52_N_to_Fe53 == ir52ng (forward), ir53gn (reverse)
            // Rate 2 == Fe53_N_to_Fe54 == ir53ng (forward), ir54gn (reverse)

            Real denom   = rates1.rr + state.y(N) * rates2.fr;
            Real denomdt = rates1.rrdt + state.y(N) * rates2.frdt;

            if (denom > 1.0e-50_rt && state.tf.t9 > 1.5_rt) {
                Real zz = 1.0_rt / denom;

                // ir2f54 in the original aprox21
                rates.fr = rates1.fr * rates2.fr * zz;
                rates.frdt = rates1.frdt * rates2.fr * zz + rates1.fr * rates2.frdt * zz - rates.fr * zz * denomdt;

                // ir1f54 in the original aprox21
                rates.rr = rates2.rr * rates1.rr * zz;
                rates.rrdt = rates2.rrdt * rates1.rr * zz + rates2.rr * rates1.rrdt * zz - rates.rr * zz * denomdt;
            }
        }

        if constexpr (rate == Fe54_2P_to_Ni56)
        {
            // In the original aprox21:
            // Rate 1 == Fe54_P_to_Co55 == irfepg (forward), ircogp (reverse)
            // Rate 2 == Co55_P_to_Ni56 == ircopg (forward), irnigp (reverse)
            // Rate 3 == Fe52_He4_to_Co55_P == irfeap (forward), ircopa (reverse)

            Real denom = rates1.rr + state.y(P) * (rates2.fr + rates3.rr);

            if (denom > 1.0e-50_rt && state.tf.t9 > 1.5_rt) {
                Real denomdt = rates1.rrdt + state.y(P) * (rates2.frdt + rates3.rrdt);

                Real zz = 1.0_rt / denom;

                // ir3f54 in the original aprox21
                rates.fr = rates1.fr * rates2.fr * zz;
                rates.frdt = rates1.frdt * rates2.fr * zz + rates1.fr * rates2.frdt * zz - rates.fr * zz * denomdt;

                // ir4f54 in the original aprox21
                rates.rr = rates2.rr * rates1.rr * zz;
                rates.rrdt = rates2.rrdt * rates1.rr * zz + rates2.rr * rates1.rrdt * zz - rates.rr * zz * denomdt;
            }
        }

        if constexpr (rate == Fe52_He4_to_Fe54_2P)
        {
            // In the original aprox21:
            // Rate 1 == Fe52_He4_to_Co55_P == irfeap (forward), ircopa (reverse)
            // Rate 2 == Fe54_P_to_Co55 == irfepg (forward), ircogp (reverse)
            // Rate 3 == Co55_P_to_Ni56 == ircopg (forward), irnigp (reverse)

            Real denom = rates2.rr + state.y(P) * (rates3.fr + rates1.rr);

            if (denom > 1.0e-50_rt && state.tf.t9 > 1.5_rt) {
                Real denomdt = rates2.rrdt + state.y(P) * (rates3.frdt + rates1.rrdt);

                Real zz = 1.0_rt / denom;

                // ir6f54 in the original aprox21
                rates.fr = rates1.fr * rates2.rr * zz;
                rates.frdt = rates1.frdt * rates2.rr * zz + rates1.fr * rates2.rrdt * zz - rates.fr * zz * denomdt;

                // ir5f54 in the original aprox21
                rates.rr = rates2.fr * rates1.rr * zz;
                rates.rrdt = rates2.frdt * rates1.rr * zz + rates2.fr * rates1.rrdt * zz - rates.rr * zz * denomdt;
            }
        }

        if constexpr (rate == Fe52_He4_to_Ni56)
        {
            // In the original aprox21:
            // Rate 1 == Fe52_He4_to_Co55_P == irfeap (forward), ircopa (reverse)
            // Rate 2 == Co55_P_to_Ni56 == ircopg (forward), irnigp (reverse)
            // Rate 3 == Fe54_P_to_Co55 == irfepg (forward), ircogp (reverse)

            Real denom = rates3.rr + state.y(P) * (rates2.fr + rates1.rr);

            if (denom > 1.0e-50_rt && state.tf.t9 > 1.5_rt) {
                Real denomdt = rates3.rrdt + state.y(P) * (rates2.frdt + rates1.rrdt);

                Real zz = 1.0_rt / denom;

                // Note that in these reactions we're adding to fr and rr
                // since this reaction already has a contribution from the
                // direct alpha capture sequence Fe52(a,g)Ni56.

                // ir7f54 in the original aprox19
                Real fr_add = rates1.fr * rates2.fr * zz;
                Real frdt_add = rates1.frdt * rates2.fr * zz + rates1.fr * rates2.frdt * zz - fr_add * zz * denomdt;
                rates.fr += state.y(P) * fr_add;
                rates.frdt += state.y(P) * frdt_add;

                // ir8f54 in the original aprox19
                Real rr_add = rates2.rr * rates1.rr * zz;
                Real rrdt_add = rates2.rrdt * rates1.rr * zz + rates2.rr * rates1.rrdt * zz - rr_add * zz * denomdt;
                rates.rr += state.y(P) * rr_add;
                rates.rrdt += state.y(P) * rrdt_add;
            }
        }

        if constexpr (rate == Fe56_to_Cr56)
        {
            // Rate 1 == Ni56_to_Fe56
            rates.fr = 1.0e-4 * rates1.fr;
        }

        if constexpr (rate == Fe54_2N_to_Fe56)
        {
            // In the original aprox21:
            // Rate 1 == Fe54_N_to_Fe55 == ir54ng (forward), ir55gn (reverse)
            // Rate 2 == Fe55_N_to_Fe56 == ir55ng (forward), ir56gn (reverse)

            Real denom   = rates1.rr + state.y(N) * rates2.fr;
            Real denomdt = rates1.rrdt + state.y(N) * rates2.frdt;

            if (denom > 1.0e-50_rt && state.tf.t9 > 1.5_rt) {
                Real zz = 1.0_rt / denom;

                // irfe56_aux2 in the original aprox21
                rates.fr = rates1.fr * rates2.fr * zz;
                rates.frdt = rates1.frdt * rates2.fr * zz + rates1.fr * rates2.frdt * zz - rates.fr * zz * denomdt;

                // irfe56_aux1 in the original aprox21
                rates.rr = rates2.rr * rates1.rr * zz;
                rates.rrdt = rates2.rrdt * rates1.rr * zz + rates2.rr * rates1.rrdt * zz - rates.rr * zz * denomdt;
            }
        }

        if constexpr (rate == Fe54_He4_to_Fe56_2P)
        {
            // In the original aprox21:
            // Rate 1 == Fe54_He4_to_Co57_P == ir54ap (forward), irco57pa (reverse)
            // Rate 2 == Fe56_P_to_Co57 == irfe56pg (forward), irco57gp (reverse)

            Real denom   = rates2.rr + state.y(P) * rates1.rr;
            Real denomdt = rates2.rrdt + state.y(P) * rates1.rrdt;

            if (denom > 1.0e-50_rt && state.tf.t9 > 1.5_rt) {
                Real zz = 1.0_rt / denom;

                // irfe56_aux4 in the original aprox21
                rates.fr = rates1.fr * rates2.rr * zz;
                rates.frdt = rates1.frdt * rates2.rr * zz + rates1.fr * rates2.rrdt * zz - rates.fr * zz * denomdt;

                // irfe56_aux3 in the original aprox21
                rates.rr = rates2.fr * rates1.rr * zz;
                rates.rrdt = rates2.frdt * rates1.rr * zz + rates2.fr * rates1.rrdt * zz - rates.rr * zz * denomdt;
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
