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

namespace table
{

  constexpr int npts = 46221;

  constexpr int ntemp = 71;
  constexpr int nden = 31;
  constexpr int nye = 21;

  extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, npts> ttlog;
  extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, npts> ddlog;
  extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, npts> yetab;

  extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, npts> abartab;
  extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, npts> ebtab;
  extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, npts> wratetab;

  extern AMREX_GPU_MANAGED amrex::Array2D<amrex::Real, 1, NumSpec, 1, npts> massfractab;

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
        Ni56_to_Fe54,
        He3_He3_to_He4_2P,
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
        NumRates=Fe54_P_to_Co55
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

            data.apply_identical_particle_factor = 0;

            data.additional_reaction_1 = Fe54_P_to_Co55;
            data.additional_reaction_2 = Co55_P_to_Ni56;
            data.additional_reaction_3 = Fe52_He4_to_Co55_P;
            break;

        case Ni56_to_Fe54:
            data.species_A = Ni56;
            data.species_D = Fe54;

            data.number_A = 54;
            data.number_D = 56;

            data.exponent_A = 1;
            data.exponent_D = 1;

            data.rate_can_be_tabulated = 0;

            data.screen_forward_reaction = 0;
            data.screen_reverse_reaction = 0;
            break;
        }

        return data;
    }

    template<int rate>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void evaluate_analytical_rate (const Real& rho, const tf_t& tf, const Real& y_e,
                                   const Real& eta, const Array1D<Real, 1, NumSpec>& y,
                                   Real& fr, Real& frdt, Real& rr, Real& rrdt)
    {
        using namespace Species;
        using namespace Rates;

        if constexpr (rate == H1_H1_to_He3) {
            Real tmp1, tmp2;
            rate_pp(tf, 1.0_rt, fr, frdt, tmp1, tmp2);
        }
        else if constexpr (rate == P_to_N) {
            if (tf.temp >= 1.0e6_rt && rho >= 1.0e-9_rt) {
                Real tmp1, tmp2;
                ecapnuc(eta, tf.temp, fr, rr, tmp1, tmp2);
            }
        }
        else if constexpr (rate == P_N_to_H2) {
            rate_png(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == H2_P_to_He3) {
            rate_dpg(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == He3_N_to_He4) {
            rate_he3ng(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == He3_He3_to_He4_2H1) {
            Real tmp1, tmp2;
            rate_he3he3(tf, 1.0_rt, fr, frdt, tmp1, tmp2);
        }
        else if constexpr (rate == He3_He4_H1_to_2He4) {
            Real tmp1, tmp2;
            rate_he3he4(tf, 1.0_rt, fr, frdt, tmp1, tmp2);
        }
        else if constexpr (rate == C12_He4_to_O16) {
            rate_c12ag(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == He4_He4_He4_to_C12) {
            rate_triplealf(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == C12_H1_to_N13) {
            rate_c12pg(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == C12_C12_to_Ne20_He4) {
            rate_c12c12(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == C12_O16_to_Mg24_He4) {
            rate_c12o16(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == N14_H1_to_O15) {
            rate_n14pg(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == N15_H1_to_O16) {
            rate_n15pg(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == N15_H1_to_C12_He4) {
            rate_n15pa(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == N14_He4_to_F18) {
            rate_n14ag(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == O16_H1_to_F17) {
            rate_o16pg(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == O16_O16_to_P31_P) {
            rate_o16o16(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == O16_He4_to_Ne20) {
            rate_o16ag(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == Ne20_He4_to_Mg24) {
            rate_ne20ag(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == Mg24_He4_to_Si28) {
            rate_mg24ag(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == Si28_He4_to_S32) {
            rate_si28ag(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == S32_He4_to_Ar36) {
            rate_s32ag(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == Ar36_He4_to_Ca40) {
            rate_ar36ag(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == Ca40_He4_to_Ti44) {
            rate_ca40ag(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == Ti44_He4_to_Cr48) {
            rate_ti44ag(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == Cr48_He4_to_Fe52) {
            rate_cr48ag(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == Fe52_He4_to_Ni56) {
            rate_fe52ag(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == Mg24_He4_to_Al27_P) {
            rate_mg24ap(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == Al27_P_to_Si28) {
            rate_al27pg(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == Si28_He4_to_P31_P) {
            rate_si28ap(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == P31_P_to_S32) {
            rate_p31pg(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == S32_He4_to_Cl35_P) {
            rate_s32ap(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == Cl35_P_to_Ar36) {
            rate_cl35pg(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == Ar36_He4_to_K39_P) {
            rate_ar36ap(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == K39_P_to_Ca40) {
            rate_k39pg(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == Ca40_He4_to_Sc43_P) {
            rate_ca40ap(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == Sc43_P_to_Ti44) {
            rate_sc43pg(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == Ti44_He4_to_V47_P) {
            rate_ti44ap(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == V47_P_to_Cr48) {
            rate_v47pg(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == Cr48_He4_to_Mn51_P) {
            rate_cr48ap(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == Mn51_P_to_Fe52) {
            rate_mn51pg(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == Fe52_He4_to_Co55_P) {
            rate_fe52ap(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == Co55_P_to_Ni56) {
            rate_co55pg(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == Fe52_N_to_Fe53) {
            rate_fe52ng(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == Fe53_N_to_Fe54) {
            rate_fe53ng(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == Fe54_P_to_Co55) {
            rate_fe54pg(tf, 1.0_rt, fr, frdt, rr, rrdt);
        }
        else if constexpr (rate == Ni56_to_Fe54) {
            if (tf.temp >= 1.0e6_rt && rho >= 1.0e-9_rt) {
                Real tmp;
                langanke(tf.temp, rho, y(Ni56), y_e, fr, tmp);

                // We scaled the number of nucleons per reaction
                // in defining the reaction so that the number of reactants
                // on each side would be integers while still conserving mass,
                // so to compensate, we scale the rate by 1/54. This serves
                // the equivalent role to the c54 constant in the original aprox21.

                fr *= 1.0_rt / 54.0_rt;
            }
        }
    }

    template<int rate>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void postprocess_rate (const Real& rho, const tf_t& tf, const Real& y_e,
                           const Real& eta, const Array1D<Real, 1, NumSpec>& y,
                           Real& fr, Real& frdt, Real& rr, Real& rrdt,
                           Real& fr1, Real& fr1dt, Real& rr1, Real& rr1dt,
                           Real& fr2, Real& fr2dt, Real& rr2, Real& rr2dt,
                           Real& fr3, Real& fr3dt, Real& rr3, Real& rr3dt)
    {
        using namespace Species;
        using namespace Rates;

        if constexpr (rate == P_P_N_N_to_He4)
        {
            // In the original aprox21:
            // Rate 1 == P_N_to_H2 == irhng (forward), irdgn (reverse)
            // Rate 2 == H2_P_to_He3 == irdpg (forward), irhegp (reverse)
            // Rate 3 == He3_N_to_He4 == irheng (forward), irhegn (reverse)

            Real denom = rr2 * rr1 + y(N) * fr3 * rr1 + y(N) * y(P) * fr3 * fr2;

            if (denom > 1.0e-50_rt && tf.t9 > 1.5_rt) {
                Real denomdt = rr2dt * rr1 + rr2 * rr1dt +
                               y(N) * (fr3dt * rr1 + fr3 * rr1dt) +
                               y(N) * y(P) * (fr3dt * fr2 + fr3 * fr2dt);

                Real zz = 1.0_rt / denom;

                // iralf1 in the original aprox21
                rr = rr3 * rr2 * rr1 * zz;
                rrdt = rr3dt * rr2 * rr1 * zz +
                       rr3 * rr2dt * rr1 * zz +
                       rr3 * rr2 * rr1dt * zz -
                       rr * zz * denomdt;

                // iralf2 in the original aprox21
                fr = fr3 * fr2 * fr1 * zz;
                frdt = fr3dt * fr2 * fr1 * zz +
                       fr3 * fr2dt * fr1 * zz +
                       fr3 * fr2 * fr1dt * zz -
                       fr * zz * denomdt;
            }
        }

        // Set this equal to the p+e->n rate.

        if constexpr (rate == H1_H1_H1_to_He3) {
            fr = fr1;
        }

        // Beta limit He3 + He4 by the B8 decay half-life

        if constexpr (rate == He3_He4_H1_to_2He4)
        {
            Real xx = 0.896_rt / y(He4);
            fr = amrex::min(fr, xx);
            if (fr == xx) {
                frdt = 0.0_rt;
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
            Real denom      = rr1 + fr2;
            Real denomdt    = rr1dt + fr2dt;

            if (denom > 1.0e-30_rt) {
                Real zz = 1.0_rt / denom;
                rate_ir    = rr1 * zz;
                dratedt_ir = (rr1dt - rate_ir * denomdt) * zz;
            }

            fr += fr1 * (1.0_rt - rate_ir);
            frdt += fr1dt * (1.0_rt - rate_ir) - fr1 * dratedt_ir;

            rr += rr2 * rate_ir;
            rrdt += rr2dt * rate_ir + rr2 * dratedt_ir;
        }

        // The O16+O16->Si28+He4 reaction has an additional contribution from the O16+O16->P31+P
        // reaction followed by the reverse Si28+He4<->P31+P reaction.
        if constexpr (rate == O16_O16_to_Si28_He4)
        {
            Real rate_ir    = 0.0e0_rt;
            Real dratedt_ir = 0.0e0_rt;
            Real denom      = rr1 + fr2;
            Real denomdt    = rr1dt + fr2dt;

            if (denom > 1.0e-30_rt) {
                Real zz = 1.0_rt / denom;
                rate_ir    = rr1 * zz;
                dratedt_ir = (rr1dt - rate_ir * denomdt) * zz;
            }

            fr = 0.56_rt * fr3 + 0.34_rt * fr3 * rate_ir;
            frdt = 0.56_rt * fr3dt + 0.34_rt * fr3dt * rate_ir + 0.34_rt * fr3 * dratedt_ir;
        }

        // The O16+O16->S32 reaction has an additional contribution from the O16+O16->P31+P
        // reaction followed by the forward P31+P<->S32 reaction.
        if constexpr (rate == O16_O16_to_S32)
        {
            Real rate_ir    = 0.0e0_rt;
            Real dratedt_ir = 0.0e0_rt;
            Real denom      = rr1 + fr2;
            Real denomdt    = rr1dt + fr2dt;

            if (denom > 1.0e-30_rt) {
                Real zz = 1.0_rt / denom;
                rate_ir    = rr1 * zz;
                dratedt_ir = (rr1dt - rate_ir * denomdt) * zz;
            }

            fr = 0.1_rt * fr3 + 0.34_rt * fr3 * (1.0_rt - rate_ir);
            frdt = 0.1_rt * fr3dt + 0.34_rt * fr3dt * (1.0_rt - rate_ir) - 0.34_rt * fr3 * dratedt_ir;
        }

        if constexpr (rate == C12_2H1_to_N14)
        {
            fr = fr1;
            frdt = fr1dt;
        }

        if constexpr (rate == C12_O16_to_Si28)
        {
            fr = fr1;
            frdt = fr1dt;
            rr = rr1;
            rrdt = rr1dt;
        }

        if constexpr (rate == N14_2H1_to_O16)
        {
            // Rate 1: N14_H1_to_O15

            // Beta limit

            if (y(H1) > 1.0e-30_rt) {
                Real xx = 5.68e-3_rt / (y(H1) * 1.57_rt);
                fr1 = amrex::min(fr1, xx);
                if (fr1 == xx) {
                    fr1dt = 0.0_rt;
                }
            }

            // Rate 2: O15_H1_to_O16
            // Rate 3: O15_H1_to_C12_He4

            Real tot = fr2 + fr3;
            Real invtot = 1.0_rt / tot;

            fr = fr1 * (fr2 * invtot);

            Real dtotdt = fr2dt + fr3dt;

            frdt = fr1 * (fr2dt * invtot - fr2 * invtot * invtot * dtotdt);
        }

        if constexpr (rate == N14_2H1_to_C12_He4)
        {
            // Rate 1: N14_H1_to_O15

            // Beta limit

            if (y(H1) > 1.0e-30_rt) {
                Real xx = 5.68e-3_rt / (y(H1) * 1.57_rt);
                fr1 = amrex::min(fr1, xx);
                if (fr1 == xx) {
                    fr1dt = 0.0_rt;
                }
            }

            // Rate 2: O15_H1_to_O16
            // Rate 3: O15_H1_to_C12_He4

            Real tot = fr2 + fr3;
            Real invtot = 1.0_rt / tot;

            fr = fr1 * (fr3 * invtot);

            Real dtotdt = fr2dt + fr3dt;

            frdt = fr1 * (fr3dt * invtot - fr3 * invtot * invtot * dtotdt);
        }

        if constexpr (rate == N14_He4_to_Ne20)
        {
            // We doubled the number of nucleons per reaction
            // in defining the reaction so that the number of reactants
            // on each side would be integers, so to compensate,
            // we'll halve the rate.

            fr = fr1 * 0.5_rt;
            frdt = fr1dt * 0.5_rt;
        }

        if constexpr (rate == O16_2H1_to_N14_He4)
        {
            // Beta limit

            Real xx = 0.0105_rt / y(H1);
            fr1 = amrex::min(fr1, xx);
            if (fr1 == xx) {
                fr1dt = 0.0_rt;
            }

            fr = fr1;
            frdt = fr1dt;
        }

        if constexpr (rate == Fe52_2N_to_Fe54)
        {
            // In the original aprox21:
            // Rate 1 == Fe52_N_to_Fe53 == ir52ng (forward), ir53gn (reverse)
            // Rate 2 == Fe53_N_to_Fe54 == ir53ng (forward), ir54gn (reverse)

            Real denom   = rr1 + y(N) * fr2;
            Real denomdt = rr1dt + y(N) * fr2dt;

            if (denom > 1.0e-50_rt && tf.t9 > 1.5_rt) {
                Real zz = 1.0_rt / denom;

                // ir2f54 in the original aprox21
                fr = fr1 * fr2 * zz;
                frdt = fr1dt * fr2 * zz + fr1 * fr2dt * zz - fr * zz * denomdt;

                // ir1f54 in the original aprox21
                rr = rr2 * rr1 * zz;
                rrdt = rr2dt * rr1 * zz + rr2 * rr1dt * zz - rr * zz * denomdt;
            }
        }

        if constexpr (rate == Fe54_2P_to_Ni56)
        {
            // In the original aprox21:
            // Rate 1 == Fe54_P_to_Co55 == irfepg (forward), ircogp (reverse)
            // Rate 2 == Co55_P_to_Ni56 == ircopg (forward), irnigp (reverse)
            // Rate 3 == Fe52_He4_to_Co55_P == irfeap (forward), ircopa (reverse)

            Real denom = rr1 + y(P) * (fr2 + rr3);

            if (denom > 1.0e-50_rt && tf.t9 > 1.5_rt) {
                Real denomdt = rr1dt + y(P) * (fr2dt + rr3dt);

                Real zz = 1.0_rt / denom;

                // ir3f54 in the original aprox21
                fr = fr1 * fr2 * zz;
                frdt = fr1dt * fr2 * zz + fr2 * fr1dt * zz - fr * zz * denomdt;

                // ir4f54 in the original aprox21
                rr = rr2 * rr1 * zz;
                rrdt = rr2dt * rr1 * zz + rr2 * rr1dt * zz - rr * zz * denomdt;
            }
        }

        if constexpr (rate == Fe52_He4_to_Fe54_2P)
        {
            // In the original aprox21:
            // Rate 1 == Fe52_He4_to_Co55_P == irfeap (forward), ircopa (reverse)
            // Rate 2 == Fe54_P_to_Co55 == irfepg (forward), ircogp (reverse)
            // Rate 3 == Co55_P_to_Ni56 == ircopg (forward), irnigp (reverse)

            Real denom = rr2 + y(P) * (fr3 + rr1);

            if (denom > 1.0e-50_rt && tf.t9 > 1.5_rt) {
                Real denomdt = rr2dt + y(P) * (fr3dt + rr1dt);

                Real zz = 1.0_rt / denom;

                // ir6f54 in the original aprox21
                fr = fr1 * rr2 * zz;
                frdt = fr1dt * rr2 * zz + fr1 * rr2dt * zz - fr * zz * denomdt;

                // ir5f54 in the original aprox21
                rr = fr2 * rr1 * zz;
                rrdt = fr2dt * rr1 * zz + fr2 * rr1dt * zz - rr * zz * denomdt;
            }
        }

        if constexpr (rate == Fe52_He4_to_Ni56)
        {
            // In the original aprox21:
            // Rate 1 == Fe52_He4_to_Co55_P == irfeap (forward), ircopa (reverse)
            // Rate 2 == Co55_P_to_Ni56 == ircopg (forward), irnigp (reverse)
            // Rate 3 == Fe54_P_to_Co55 == irfepg (forward), ircogp (reverse)

            Real denom = rr3 + y(P) * (fr2 + rr1);

            if (denom > 1.0e-50_rt && tf.t9 > 1.5_rt) {
                Real denomdt = rr3dt + y(P) * (fr2dt + rr1dt);

                Real zz = 1.0_rt / denom;

                // Note that in these reactions we're adding to fr and rr
                // since this reaction already has a contribution from the
                // direct alpha capture sequence Fe52(a,g)Ni56.

                // ir7f54 in the original aprox19
                Real fr_add = fr1 * fr2 * zz;
                Real frdt_add = fr1dt * fr2 * zz + fr1 * fr2dt * zz - fr_add * zz * denomdt;
                fr += y(P) * fr_add;
                frdt += y(P) * frdt_add;

                // ir8f54 in the original aprox19
                Real rr_add = rr2 * rr1 * zz;
                Real rrdt_add = rr2dt * rr1 * zz + rr2 * rr1dt * zz - rr_add * zz * denomdt; 
                rr += y(P) * rr_add;
                rrdt += y(P) * rrdt_add;
            }
        }
    }

} // namespace RHS

#endif
