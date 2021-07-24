#ifndef _actual_network_H_
#define _actual_network_H_

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>

#include <fundamental_constants.H>
#include <network_properties.H>
#include <network_utilities.H>

using namespace amrex;

void actual_network_init();

const std::string network_name = "aprox13";

namespace network
{
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> bion;
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> mion;
}

namespace Rates {
    enum NetworkRates {
                        He4_He4_He4_to_C12_forward = 1,
                        He4_He4_He4_to_C12_reverse,
                        C12_He4_to_O16_forward,
                        C12_He4_to_O16_reverse,
                        C12_C12_to_Ne20_He4_forward,
                        C12_C12_to_Ne20_He4_reverse,
                        C12_O16_to_Mg24_He4_forward,
                        C12_O16_to_Mg24_He4_reverse,
                        C12_O16_to_Si28_forward,
                        C12_O16_to_Si28_reverse,
                        O16_O16_to_Si28_He4_forward,
                        O16_O16_to_Si28_He4_reverse,
                        O16_O16_to_S32_forward,
                        O16_O16_to_S32_reverse,
                        O16_He4_to_Ne20_forward,
                        O16_He4_to_Ne20_reverse,
                        Ne20_He4_to_Mg24_forward,
                        Ne20_He4_to_Mg24_reverse,
                        Mg24_He4_to_Si28_forward,
                        Mg24_He4_to_Si28_reverse,
                        Si28_He4_to_S32_forward,
                        Si28_He4_to_S32_reverse,
                        S32_He4_to_Ar36_forward,
                        S32_He4_to_Ar36_reverse,
                        Ar36_He4_to_Ca40_forward,
                        Ar36_He4_to_Ca40_reverse,
                        Ca40_He4_to_Ti44_forward,
                        Ca40_He4_to_Ti44_reverse,
                        Ti44_He4_to_Cr48_forward,
                        Ti44_He4_to_Cr48_reverse,
                        Cr48_He4_to_Fe52_forward,
                        Cr48_He4_to_Fe52_reverse,
                        Fe52_He4_to_Ni56_forward,
                        Fe52_He4_to_Ni56_reverse,
                        Mg24_He4_to_Al27_P_forward,
                        Mg24_He4_to_Al27_P_reverse,
                        Al27_P_to_Si28_forward,
                        Al27_P_to_Si28_reverse,
                        Si28_He4_to_P31_P_forward,
                        Si28_He4_to_P31_P_reverse,
                        P31_P_to_S32_forward,
                        P31_P_to_S32_reverse,
                        S32_He4_to_Cl35_P_forward,
                        S32_He4_to_Cl35_P_reverse,
                        Cl35_P_to_Ar36_forward,
                        Cl35_P_to_Ar36_reverse,
                        Ar36_He4_to_K39_P_forward,
                        Ar36_He4_to_K39_P_reverse,
                        K39_P_to_Ca40_forward,
                        K39_P_to_Ca40_reverse,
                        Ca40_He4_to_Sc43_P_forward,
                        Ca40_He4_to_Sc43_P_reverse,
                        Sc43_P_to_Ti44_forward,
                        Sc43_P_to_Ti44_reverse,
                        Ti44_He4_to_V47_P_forward,
                        Ti44_He4_to_V47_P_reverse,
                        V47_P_to_Cr48_forward,
                        V47_P_to_Cr48_reverse,
                        Cr48_He4_to_Mn51_P_forward,
                        Cr48_He4_to_Mn51_P_reverse,
                        Mn51_P_to_Fe52_forward,
                        Mn51_P_to_Fe52_reverse,
                        Fe52_He4_to_Co55_P_forward,
                        Fe52_He4_to_Co55_P_reverse,
                        Co55_P_to_Ni56_forward,
                        Co55_P_to_Ni56_reverse,
                        NumRates=Co55_P_to_Ni56_reverse
                      };

    enum NetworkRatesFR {
        He4_He4_He4_to_C12 = 1,
        C12_He4_to_O16,
        C12_C12_to_Ne20_He4,
        C12_O16_to_Mg24_He4,
        C12_O16_to_Si28,
        O16_O16_to_Si28_He4,
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
        NumRatesFR=Fe52_He4_to_Ni56
    };

    const int NumGroups = 2;
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

        case O16_O16_to_S32:
            data.species_A = O16;
            data.species_D = S32;

            data.number_A = 2;
            data.number_D = 1;

            data.exponent_A = 2;
            data.exponent_D = 1;
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
            break;

        }

        return data;
    }

}

#endif
