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

const std::string network_name = "iso7";

namespace network
{
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> bion;
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> mion;
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> wion;
}

namespace Rates {
    enum NetworkRatesFR {
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
        NumRatesFR = Si28_7He4_to_Ni56
    };

    enum NetworkRates {C12_He4_to_O16_forward = 1,
                       C12_He4_to_O16_reverse,
                       He4_He4_He4_to_C12_forward,
                       He4_He4_He4_to_C12_reverse,
                       C12_C12_to_Ne20_He4_forward,
                       C12_C12_to_Ne20_He4_reverse,
                       C12_O16_to_Mg24_He4_forward,
                       C12_O16_to_Mg24_He4_reverse,
                       C12_O16_to_Si28_forward,
                       C12_O16_to_Si28_reverse,
                       O16_O16_to_Si28_He4_forward,
                       O16_O16_to_Si28_He4_reverse,
                       O16_He4_to_Ne20_forward,
                       O16_He4_to_Ne20_reverse,
                       Ne20_He4_to_Mg24_forward,
                       Ne20_He4_to_Mg24_reverse,
                       Mg24_He4_to_Si28_forward,
                       Mg24_He4_to_Si28_reverse,
                       Ca40_He4_to_Ti44_forward,
                       Ca40_He4_to_Ti44_reverse,
                       Si28_7He4_to_Ni56_forward,
                       Si28_7He4_to_Ni56_reverse,
                       NumRates = Si28_7He4_to_Ni56_reverse};

    const int NumGroups = 2 + NumSpec;
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
            break;

        }

        return data;
    }

    // Contribution to the Jacobian from d(RHS(species1)) / d(species2)
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    constexpr jac_t jac_data (int species1, int species2, int rateindex)
    {
        using namespace Species;
        using namespace Rates;

        jac_t data{};

        // The RHS for species 1 has the form prefactor * Y(speciesA) * Y(speciesB) * Y(speciesC) * rate.
        // If species2 is one of species{A,B,C}, then the Jacobian has a contribution from the derivative
        // of the RHS term with respect to Y(species2). Additionally, the rate itself may have a functional
        // dependence on species2 so we also return the same format for d(rate) / d(species2) if it exists.

        data.prefactor  = 0.0_rt;
        data.specindex1 = -1;
        data.specindex2 = -1;
        data.specindex3 = -1;

        data.rate_prefactor = 0.0_rt;
        data.rate_specindex1 = -1;
        data.rate_specindex2 = -1;
        data.rate_specindex3 = -1;

        switch (jac_rate<NumSpec, NumRates>(species1, species2, rateindex)) {

        case jac_rate<NumSpec, NumRates>(He4, He4, Si28_7He4_to_Ni56_forward):
            data.rate_prefactor = -7.0_rt;
            data.rate_specindex1 = He4;
            data.rate_specindex2 = Si28;
            break;

        case jac_rate<NumSpec, NumRates>(He4, He4, Si28_7He4_to_Ni56_reverse):
            data.rate_prefactor = 7.0_rt;
            data.rate_specindex1 = Ni56;
            break;

        case jac_rate<NumSpec, NumRates>(Si28, He4, Si28_7He4_to_Ni56_forward):
            data.rate_prefactor = -1.0_rt;
            data.rate_specindex1 = He4;
            data.rate_specindex2 = Si28;
            break;

        case jac_rate<NumSpec, NumRates>(Si28, He4, Si28_7He4_to_Ni56_reverse):
            data.rate_prefactor = 1.0_rt;
            data.rate_specindex1 = Ni56;
            break;

        case jac_rate<NumSpec, NumRates>(Ni56, He4, Si28_7He4_to_Ni56_forward):
            data.rate_prefactor = 1.0_rt;
            data.rate_specindex1 = He4;
            data.rate_specindex2 = Si28;
            break;

        case jac_rate<NumSpec, NumRates>(Ni56, He4, Si28_7He4_to_Ni56_reverse):
            data.rate_prefactor = -1.0_rt;
            data.rate_specindex1 = Ni56;
            break;

        }

        return data;
    }
}

#endif
