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

namespace C
{
    namespace Legacy
    {
        // These are the values of the constants used in the original iso7
        constexpr amrex::Real m_n = 1.67492721184e-24_rt;
        constexpr amrex::Real m_p = 1.67262163783e-24_rt;
        constexpr amrex::Real m_e = 9.1093821545e-28_rt;

        constexpr amrex::Real eV2erg  = 1.60217648740e-12_rt;
        constexpr amrex::Real MeV2erg = eV2erg*1.0e6_rt;
        constexpr amrex::Real MeV2gr  = MeV2erg/(c_light*c_light);

        constexpr amrex::Real n_A = 6.0221417930e23_rt;

        // conversion factor for nuclear energy generation rate
        constexpr amrex::Real enuc_conv2 = -n_A * c_light * c_light;
    }
}

const std::string network_name = "iso7";

namespace iso7
{
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> bion;
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> mion;
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> wion;
}

namespace Rates {
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
    constexpr rhs_t rhs_data (int species, int rateindex)
    {
        using namespace Species;
        using namespace Rates;

        rhs_t data{};

        data.prefactor  = 0.0_rt;
        data.specindex1 = -1;
        data.specindex2 = -1;
        data.specindex3 = -1;

        switch (rhs_rate<NumRates>(species, rateindex)) {

        case rhs_rate<NumRates>(He4, He4_He4_He4_to_C12_reverse):
            data.prefactor  = 3.0_rt;
            data.specindex1 = C12;
            break;

        case rhs_rate<NumRates>(He4, He4_He4_He4_to_C12_forward):
            data.prefactor  = -3.0_rt;
            data.specindex1 = He4;
            data.specindex2 = He4;
            data.specindex3 = He4;
            break;

        case rhs_rate<NumRates>(He4, C12_He4_to_O16_reverse):
            data.prefactor  = 1.0_rt;
            data.specindex1 = O16;
            break;

        case rhs_rate<NumRates>(He4, C12_He4_to_O16_forward):
            data.prefactor  = -1.0_rt;
            data.specindex1 = C12;
            data.specindex2 = He4;
            break;

        case rhs_rate<NumRates>(He4, C12_C12_to_Ne20_He4_forward):
            data.prefactor  = 1.0_rt;
            data.specindex1 = C12;
            data.specindex2 = C12;
            break;

        case rhs_rate<NumRates>(He4, C12_O16_to_Mg24_He4_forward):
            data.prefactor  = 1.0_rt;
            data.specindex1 = C12;
            data.specindex2 = O16;
            break;

        case rhs_rate<NumRates>(He4, O16_O16_to_Si28_He4_forward):
            data.prefactor  = 1.0_rt;
            data.specindex1 = O16;
            data.specindex2 = O16;
            break;

        case rhs_rate<NumRates>(He4, O16_He4_to_Ne20_forward):
            data.prefactor  = -1.0_rt;
            data.specindex1 = O16;
            data.specindex2 = He4;
            break;

        case rhs_rate<NumRates>(He4, O16_He4_to_Ne20_reverse):
            data.prefactor  = 1.0_rt;
            data.specindex1 = Ne20;
            break;

        case rhs_rate<NumRates>(He4, Ne20_He4_to_Mg24_reverse):
            data.prefactor  = 1.0_rt;
            data.specindex1 = Mg24;
            break;

        case rhs_rate<NumRates>(He4, Ne20_He4_to_Mg24_forward):
            data.prefactor  = -1.0_rt;
            data.specindex1 = Ne20;
            data.specindex2 = He4;
            break;

        case rhs_rate<NumRates>(He4, Mg24_He4_to_Si28_reverse):
            data.prefactor  = 1.0_rt;
            data.specindex1 = Si28;
            break;

        case rhs_rate<NumRates>(He4, Mg24_He4_to_Si28_forward):
            data.prefactor  = -1.0_rt;
            data.specindex1 = Mg24;
            data.specindex2 = He4;
            break;

        case rhs_rate<NumRates>(He4, Si28_7He4_to_Ni56_forward):
            data.prefactor  = -7.0_rt;
            data.specindex1 = He4;
            data.specindex2 = Si28;
            break;

        case rhs_rate<NumRates>(He4, Si28_7He4_to_Ni56_reverse):
            data.prefactor  = 7.0_rt;
            data.specindex1 = Ni56;
            break;

        case rhs_rate<NumRates>(C12, He4_He4_He4_to_C12_forward):
            data.prefactor  = 1.0_rt;
            data.specindex1 = He4;
            data.specindex2 = He4;
            data.specindex3 = He4;
            break;

        case rhs_rate<NumRates>(C12, He4_He4_He4_to_C12_reverse):
            data.prefactor  = -1.0_rt;
            data.specindex1 = C12;
            break;

        case rhs_rate<NumRates>(C12, C12_He4_to_O16_reverse):
            data.prefactor  = 1.0_rt;
            data.specindex1 = O16;
            break;

        case rhs_rate<NumRates>(C12, C12_He4_to_O16_forward):
            data.prefactor  = -1.0_rt;
            data.specindex1 = C12;
            data.specindex2 = He4;
            break;

        case rhs_rate<NumRates>(C12, C12_C12_to_Ne20_He4_forward):
            data.prefactor  = -2.0_rt;
            data.specindex1 = C12;
            data.specindex2 = C12;
            break;

        case rhs_rate<NumRates>(C12, C12_O16_to_Mg24_He4_forward):
            data.prefactor  = -1.0_rt;
            data.specindex1 = C12;
            data.specindex2 = O16;
            break;

        case rhs_rate<NumRates>(C12, C12_O16_to_Si28_forward):
            data.prefactor  = -1.0_rt;
            data.specindex1 = C12;
            data.specindex2 = O16;
            break;

        case rhs_rate<NumRates>(O16, C12_He4_to_O16_reverse):
            data.prefactor  = -1.0_rt;
            data.specindex1 = O16;
            break;

        case rhs_rate<NumRates>(O16, C12_He4_to_O16_forward):
            data.prefactor  = 1.0_rt;
            data.specindex1 = C12;
            data.specindex2 = He4;
            break;

        case rhs_rate<NumRates>(O16, C12_O16_to_Mg24_He4_forward):
            data.prefactor  = -1.0_rt;
            data.specindex1 = C12;
            data.specindex2 = O16;
            break;

        case rhs_rate<NumRates>(O16, C12_O16_to_Si28_forward):
            data.prefactor  = -1.0_rt;
            data.specindex1 = C12;
            data.specindex2 = O16;
            break;

        case rhs_rate<NumRates>(O16, O16_O16_to_Si28_He4_forward):
            data.prefactor  = -2.0_rt;
            data.specindex1 = O16;
            data.specindex2 = O16;
            break;

        case rhs_rate<NumRates>(O16, O16_He4_to_Ne20_forward):
            data.prefactor  = -1.0_rt;
            data.specindex1 = O16;
            data.specindex2 = He4;
            break;

        case rhs_rate<NumRates>(O16, O16_He4_to_Ne20_reverse):
            data.prefactor  = 1.0_rt;
            data.specindex1 = Ne20;
            break;

        case rhs_rate<NumRates>(Ne20, C12_C12_to_Ne20_He4_forward):
            data.prefactor  = 1.0_rt;
            data.specindex1 = C12;
            data.specindex2 = C12;
            break;

        case rhs_rate<NumRates>(Ne20, O16_He4_to_Ne20_forward):
            data.prefactor  = 1.0_rt;
            data.specindex1 = O16;
            data.specindex2 = He4;
            break;

        case rhs_rate<NumRates>(Ne20, O16_He4_to_Ne20_reverse):
            data.prefactor  = -1.0_rt;
            data.specindex1 = Ne20;
            break;

        case rhs_rate<NumRates>(Ne20, Ne20_He4_to_Mg24_reverse):
            data.prefactor  = 1.0_rt;
            data.specindex1 = Mg24;
            break;

        case rhs_rate<NumRates>(Ne20, Ne20_He4_to_Mg24_forward):
            data.prefactor  = -1.0_rt;
            data.specindex1 = Ne20;
            data.specindex2 = He4;
            break;

        case rhs_rate<NumRates>(Mg24, C12_O16_to_Mg24_He4_forward):
            data.prefactor  = 1.0_rt;
            data.specindex1 = C12;
            data.specindex2 = O16;
            break;

        case rhs_rate<NumRates>(Mg24, Ne20_He4_to_Mg24_reverse):
            data.prefactor  = -1.0_rt;
            data.specindex1 = Mg24;
            break;

        case rhs_rate<NumRates>(Mg24, Ne20_He4_to_Mg24_forward):
            data.prefactor  = 1.0_rt;
            data.specindex1 = Ne20;
            data.specindex2 = He4;
            break;

        case rhs_rate<NumRates>(Mg24, Mg24_He4_to_Si28_reverse):
            data.prefactor  = 1.0_rt;
            data.specindex1 = Si28;
            break;

        case rhs_rate<NumRates>(Mg24, Mg24_He4_to_Si28_forward):
            data.prefactor  = -1.0_rt;
            data.specindex1 = Mg24;
            data.specindex2 = He4;
            break;

        case rhs_rate<NumRates>(Si28, C12_O16_to_Si28_forward):
            data.prefactor  = 1.0_rt;
            data.specindex1 = C12;
            data.specindex2 = O16;
            break;

        case rhs_rate<NumRates>(Si28, O16_O16_to_Si28_He4_forward):
            data.prefactor  = 1.0_rt;
            data.specindex1 = O16;
            data.specindex2 = O16;
            break;

        case rhs_rate<NumRates>(Si28, Mg24_He4_to_Si28_reverse):
            data.prefactor  = -1.0_rt;
            data.specindex1 = Si28;
            break;

        case rhs_rate<NumRates>(Si28, Mg24_He4_to_Si28_forward):
            data.prefactor  = 1.0_rt;
            data.specindex1 = Mg24;
            data.specindex2 = He4;
            break;

        case rhs_rate<NumRates>(Si28, Si28_7He4_to_Ni56_forward):
            data.prefactor  = -1.0_rt;
            data.specindex1 = He4;
            data.specindex2 = Si28;
            break;

        case rhs_rate<NumRates>(Si28, Si28_7He4_to_Ni56_reverse):
            data.prefactor  = 1.0_rt;
            data.specindex1 = Ni56;
            break;

        case rhs_rate<NumRates>(Ni56, Si28_7He4_to_Ni56_forward):
            data.prefactor  = 1.0_rt;
            data.specindex1 = He4;
            data.specindex2 = Si28;
            break;

        case rhs_rate<NumRates>(Ni56, Si28_7He4_to_Ni56_reverse):
            data.prefactor  = -1.0_rt;
            data.specindex1 = Ni56;
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

        case jac_rate<NumSpec, NumRates>(He4, He4, He4_He4_He4_to_C12_forward):
            data.prefactor = -9.0_rt;
            data.specindex1 = He4;
            data.specindex2 = He4;
            break;

        case jac_rate<NumSpec, NumRates>(He4, He4, C12_He4_to_O16_forward):
            data.prefactor = -1.0_rt;
            data.specindex1 = C12;
            break;

        case jac_rate<NumSpec, NumRates>(He4, He4, O16_He4_to_Ne20_forward):
            data.prefactor = -1.0_rt;
            data.specindex1 = O16;
            break;

        case jac_rate<NumSpec, NumRates>(He4, He4, Ne20_He4_to_Mg24_forward):
            data.prefactor = -1.0_rt;
            data.specindex1 = Ne20;
            break;

        case jac_rate<NumSpec, NumRates>(He4, He4, Mg24_He4_to_Si28_forward):
            data.prefactor = -1.0_rt;
            data.specindex1 = Mg24;
            break;

        case jac_rate<NumSpec, NumRates>(He4, He4, Si28_7He4_to_Ni56_forward):
            data.prefactor = -7.0_rt;
            data.specindex1 = Si28;

            data.rate_prefactor = -7.0_rt;
            data.rate_specindex1 = He4;
            data.rate_specindex2 = Si28;
            break;

        case jac_rate<NumSpec, NumRates>(He4, He4, Si28_7He4_to_Ni56_reverse):
            data.rate_prefactor = 7.0_rt;
            data.rate_specindex1 = Ni56;
            break;

        case jac_rate<NumSpec, NumRates>(He4, C12, He4_He4_He4_to_C12_reverse):
            data.prefactor = 3.0_rt;
            break;

        case jac_rate<NumSpec, NumRates>(He4, C12, C12_He4_to_O16_forward):
            data.prefactor = -1.0_rt;
            data.specindex1 = He4;
            break;

        case jac_rate<NumSpec, NumRates>(He4, C12, C12_C12_to_Ne20_He4_forward):
            data.prefactor = 2.0_rt;
            data.specindex1 = C12;
            break;

        case jac_rate<NumSpec, NumRates>(He4, C12, C12_O16_to_Mg24_He4_forward):
            data.prefactor = 1.0_rt;
            data.specindex1 = O16;
            break;

        case jac_rate<NumSpec, NumRates>(He4, O16, C12_He4_to_O16_reverse):
            data.prefactor = 1.0_rt;
            break;

        case jac_rate<NumSpec, NumRates>(He4, O16, C12_O16_to_Mg24_He4_forward):
            data.prefactor = 1.0_rt;
            data.specindex1 = C12;
            break;

        case jac_rate<NumSpec, NumRates>(He4, O16, O16_O16_to_Si28_He4_forward):
            data.prefactor = 2.0_rt;
            data.specindex1 = O16;
            break;

        case jac_rate<NumSpec, NumRates>(He4, O16, O16_He4_to_Ne20_forward):
            data.prefactor = -1.0_rt;
            data.specindex1 = He4;
            break;

        case jac_rate<NumSpec, NumRates>(He4, Ne20, O16_He4_to_Ne20_reverse):
            data.prefactor = 1.0_rt;
            break;

        case jac_rate<NumSpec, NumRates>(He4, Ne20, Ne20_He4_to_Mg24_forward):
            data.prefactor = -1.0_rt;
            data.specindex1 = He4;
            break;

        case jac_rate<NumSpec, NumRates>(He4, Mg24, Ne20_He4_to_Mg24_reverse):
            data.prefactor = 1.0_rt;
            break;

        case jac_rate<NumSpec, NumRates>(He4, Mg24, Mg24_He4_to_Si28_forward):
            data.prefactor = -1.0_rt;
            data.specindex1 = He4;
            break;

        case jac_rate<NumSpec, NumRates>(He4, Si28, Mg24_He4_to_Si28_reverse):
            data.prefactor = 1.0_rt;
            break;

        case jac_rate<NumSpec, NumRates>(He4, Si28, Si28_7He4_to_Ni56_forward):
            data.prefactor = -7.0_rt;
            data.specindex1 = He4;
            break;

        case jac_rate<NumSpec, NumRates>(He4, Ni56, Si28_7He4_to_Ni56_reverse):
            data.prefactor = 7.0_rt;
            break;

        case jac_rate<NumSpec, NumRates>(C12, He4, He4_He4_He4_to_C12_forward):
            data.prefactor = 3.0_rt;
            data.specindex1 = He4;
            data.specindex2 = He4;
            break;

        case jac_rate<NumSpec, NumRates>(C12, He4, C12_He4_to_O16_forward):
            data.prefactor = -1.0_rt;
            data.specindex1 = C12;
            break;

        case jac_rate<NumSpec, NumRates>(C12, C12, He4_He4_He4_to_C12_reverse):
            data.prefactor = -1.0_rt;
            break;

        case jac_rate<NumSpec, NumRates>(C12, C12, C12_He4_to_O16_forward):
            data.prefactor = -1.0_rt;
            data.specindex1 = He4;
            break;

        case jac_rate<NumSpec, NumRates>(C12, C12, C12_C12_to_Ne20_He4_forward):
            data.prefactor = -4.0_rt;
            data.specindex1 = C12;
            break;

        case jac_rate<NumSpec, NumRates>(C12, C12, C12_O16_to_Mg24_He4_forward):
            data.prefactor = -1.0_rt;
            data.specindex1 = O16;
            break;

        case jac_rate<NumSpec, NumRates>(C12, C12, C12_O16_to_Si28_forward):
            data.prefactor = -1.0_rt;
            data.specindex1 = O16;
            break;

        case jac_rate<NumSpec, NumRates>(C12, O16, C12_He4_to_O16_reverse):
            data.prefactor = 1.0_rt;
            break;

        case jac_rate<NumSpec, NumRates>(C12, O16, C12_O16_to_Mg24_He4_forward):
            data.prefactor = -1.0_rt;
            data.specindex1 = C12;
            break;

        case jac_rate<NumSpec, NumRates>(C12, O16, C12_O16_to_Si28_forward):
            data.prefactor = -1.0_rt;
            data.specindex1 = C12;
            break;

        case jac_rate<NumSpec, NumRates>(O16, He4, C12_He4_to_O16_forward):
            data.prefactor = 1.0_rt;
            data.specindex1 = C12;
            break;

        case jac_rate<NumSpec, NumRates>(O16, He4, O16_He4_to_Ne20_forward):
            data.prefactor = -1.0_rt;
            data.specindex1 = O16;
            break;

        case jac_rate<NumSpec, NumRates>(O16, C12, C12_He4_to_O16_forward):
            data.prefactor = 1.0_rt;
            data.specindex1 = He4;
            break;

        case jac_rate<NumSpec, NumRates>(O16, C12, C12_O16_to_Mg24_He4_forward):
            data.prefactor = -1.0_rt;
            data.specindex1 = O16;
            break;

        case jac_rate<NumSpec, NumRates>(O16, C12, C12_O16_to_Si28_forward):
            data.prefactor = -1.0_rt;
            data.specindex1 = O16;
            break;

        case jac_rate<NumSpec, NumRates>(O16, O16, C12_He4_to_O16_reverse):
            data.prefactor = -1.0_rt;
            break;

        case jac_rate<NumSpec, NumRates>(O16, O16, C12_O16_to_Mg24_He4_forward):
            data.prefactor = -1.0_rt;
            data.specindex1 = C12;
            break;

        case jac_rate<NumSpec, NumRates>(O16, O16, C12_O16_to_Si28_forward):
            data.prefactor = -1.0_rt;
            data.specindex1 = C12;
            break;

        case jac_rate<NumSpec, NumRates>(O16, O16, O16_O16_to_Si28_He4_forward):
            data.prefactor = -4.0_rt;
            data.specindex1 = O16;
            break;

        case jac_rate<NumSpec, NumRates>(O16, O16, O16_He4_to_Ne20_forward):
            data.prefactor = -1.0_rt;
            data.specindex1 = He4;
            break;

        case jac_rate<NumSpec, NumRates>(O16, Ne20, O16_He4_to_Ne20_reverse):
            data.prefactor = 1.0_rt;
            break;

        case jac_rate<NumSpec, NumRates>(Ne20, He4, O16_He4_to_Ne20_forward):
            data.prefactor = 1.0_rt;
            data.specindex1 = O16;
            break;

        case jac_rate<NumSpec, NumRates>(Ne20, He4, Ne20_He4_to_Mg24_forward):
            data.prefactor = -1.0_rt;
            data.specindex1 = Ne20;
            break;

        case jac_rate<NumSpec, NumRates>(Ne20, C12, C12_C12_to_Ne20_He4_forward):
            data.prefactor = 2.0_rt;
            data.specindex1 = C12;
            break;

        case jac_rate<NumSpec, NumRates>(Ne20, O16, O16_He4_to_Ne20_forward):
            data.prefactor = 1.0_rt;
            data.specindex1 = He4;
            break;

        case jac_rate<NumSpec, NumRates>(Ne20, Ne20, O16_He4_to_Ne20_reverse):
            data.prefactor = -1.0_rt;
            break;

        case jac_rate<NumSpec, NumRates>(Ne20, Ne20, Ne20_He4_to_Mg24_forward):
            data.prefactor = -1.0_rt;
            data.specindex1 = He4;
            break;

        case jac_rate<NumSpec, NumRates>(Ne20, Mg24, Ne20_He4_to_Mg24_reverse):
            data.prefactor = 1.0_rt;
            break;

        case jac_rate<NumSpec, NumRates>(Mg24, He4, Ne20_He4_to_Mg24_forward):
            data.prefactor = 1.0_rt;
            data.specindex1 = Ne20;
            break;

        case jac_rate<NumSpec, NumRates>(Mg24, He4, Mg24_He4_to_Si28_forward):
            data.prefactor = -1.0_rt;
            data.specindex1 = Mg24;
            break;

        case jac_rate<NumSpec, NumRates>(Mg24, C12, C12_O16_to_Mg24_He4_forward):
            data.prefactor = 1.0_rt;
            data.specindex1 = O16;
            break;

        case jac_rate<NumSpec, NumRates>(Mg24, O16, C12_O16_to_Mg24_He4_forward):
            data.prefactor = 1.0_rt;
            data.specindex1 = C12;
            break;

        case jac_rate<NumSpec, NumRates>(Mg24, Ne20, Ne20_He4_to_Mg24_forward):
            data.prefactor = 1.0_rt;
            data.specindex1 = He4;
            break;

        case jac_rate<NumSpec, NumRates>(Mg24, Mg24, Ne20_He4_to_Mg24_reverse):
            data.prefactor = -1.0_rt;
            break;

        case jac_rate<NumSpec, NumRates>(Mg24, Mg24, Mg24_He4_to_Si28_forward):
            data.prefactor = -1.0_rt;
            data.specindex1 = He4;
            break;

        case jac_rate<NumSpec, NumRates>(Mg24, Si28, Mg24_He4_to_Si28_reverse):
            data.prefactor = 1.0_rt;
            break;

        case jac_rate<NumSpec, NumRates>(Si28, He4, Mg24_He4_to_Si28_forward):
            data.prefactor = 1.0_rt;
            data.specindex1 = Mg24;
            break;

        case jac_rate<NumSpec, NumRates>(Si28, He4, Si28_7He4_to_Ni56_forward):
            data.prefactor = -1.0_rt;
            data.specindex1 = Si28;

            data.rate_prefactor = -1.0_rt;
            data.rate_specindex1 = He4;
            data.rate_specindex2 = Si28;
            break;

        case jac_rate<NumSpec, NumRates>(Si28, He4, Si28_7He4_to_Ni56_reverse):
            data.rate_prefactor = 1.0_rt;
            data.rate_specindex1 = Ni56;
            break;

        case jac_rate<NumSpec, NumRates>(Si28, C12, C12_O16_to_Si28_forward):
            data.prefactor = 1.0_rt;
            data.specindex1 = O16;
            break;

        case jac_rate<NumSpec, NumRates>(Si28, O16, O16_O16_to_Si28_He4_forward):
            data.prefactor = 2.0_rt;
            data.specindex1 = O16;
            break;

        case jac_rate<NumSpec, NumRates>(Si28, O16, C12_O16_to_Si28_forward):
            data.prefactor = 1.0_rt;
            data.specindex1 = C12;
            break;

        case jac_rate<NumSpec, NumRates>(Si28, Mg24, Mg24_He4_to_Si28_forward):
            data.prefactor = 1.0_rt;
            data.specindex1 = He4;
            break;

        case jac_rate<NumSpec, NumRates>(Si28, Si28, Mg24_He4_to_Si28_reverse):
            data.prefactor = -1.0_rt;
            break;

        case jac_rate<NumSpec, NumRates>(Si28, Si28, Si28_7He4_to_Ni56_forward):
            data.prefactor = -1.0_rt;
            data.specindex1 = He4;
            break;

        case jac_rate<NumSpec, NumRates>(Si28, Ni56, Si28_7He4_to_Ni56_reverse):
            data.prefactor = 1.0_rt;
            break;

        case jac_rate<NumSpec, NumRates>(Ni56, He4, Si28_7He4_to_Ni56_forward):
            data.prefactor = 1.0_rt;
            data.specindex1 = Si28;

            data.rate_prefactor = 1.0_rt;
            data.rate_specindex1 = He4;
            data.rate_specindex2 = Si28;
            break;

        case jac_rate<NumSpec, NumRates>(Ni56, He4, Si28_7He4_to_Ni56_reverse):
            data.rate_prefactor = -1.0_rt;
            data.rate_specindex1 = Ni56;
            break;

        case jac_rate<NumSpec, NumRates>(Ni56, Si28, Si28_7He4_to_Ni56_forward):
            data.prefactor = 1.0_rt;
            data.specindex1 = He4;
            break;

        case jac_rate<NumSpec, NumRates>(Ni56, Ni56, Si28_7He4_to_Ni56_reverse):
            data.prefactor = -1.0_rt;
            break;

        }

        return data;
    }
}

#endif
