#ifndef _actual_network_H_
#define _actual_network_H_

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>

#include <fundamental_constants.H>
#include <network_properties.H>

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
    enum NetworkRates {ircag = 1,
                       iroga,
                       ir3a,
                       irg3a,
                       ir1212,
                       ir1216,
                       ir1616,
                       iroag,
                       irnega,
                       irneag,
                       irmgga,
                       irmgag,
                       irsiga,
                       ircaag,
                       irtiga,
                       irsi2ni,
                       irni2si,
                       NumRates = irni2si};

    const int NumGroups = 4;

    extern amrex::Vector<std::string> names;
}

namespace RHS {

    struct rhs_t {
        Real prefactor;
        int specindex1;
        int specindex2;
        int specindex3;
    };

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

        switch (species) {

        case He4:

            switch (rateindex) {

            case irg3a:
                data.prefactor  = 3.0_rt;
                data.specindex1 = C12;
                break;

            case ir3a:
                data.prefactor  = -0.5_rt;
                data.specindex1 = He4;
                data.specindex2 = He4;
                data.specindex3 = He4;
                break;

            case iroga:
                data.prefactor  = 1.0_rt;
                data.specindex1 = O16;
                break;

            case ircag:
                data.prefactor  = -1.0_rt;
                data.specindex1 = C12;
                data.specindex2 = He4;
                break;

            case ir1212:
                data.prefactor  = 0.5_rt;
                data.specindex1 = C12;
                data.specindex2 = C12;
                break;

            case ir1216:
                data.prefactor  = 0.5_rt;
                data.specindex1 = C12;
                data.specindex2 = O16;
                break;

            case ir1616:
                data.prefactor  = 0.5_rt;
                data.specindex1 = O16;
                data.specindex2 = O16;
                break;

            case iroag:
                data.prefactor  = -1.0_rt;
                data.specindex1 = O16;
                data.specindex2 = He4;
                break;

            case irnega:
                data.prefactor  = 1.0_rt;
                data.specindex1 = Ne20;
                break;

            case irmgga:
                data.prefactor  = 1.0_rt;
                data.specindex1 = Mg24;
                break;

            case irneag:
                data.prefactor  = -1.0_rt;
                data.specindex1 = Ne20;
                data.specindex2 = He4;
                break;

            case irsiga:
                data.prefactor  = 1.0_rt;
                data.specindex1 = Si28;
                break;

            case irmgag:
                data.prefactor  = -1.0_rt;
                data.specindex1 = Mg24;
                data.specindex2 = He4;
                break;

            case irsi2ni:
                data.prefactor  = -7.0_rt;
                data.specindex1 = He4;
                break;

            case irni2si:
                data.prefactor  = 7.0_rt;
                data.specindex1 = Ni56;
                break;

            }
            break;

        case C12:

            switch (rateindex) {

            case ir3a:
                data.prefactor  = 1.0_rt / 6.0_rt;
                data.specindex1 = He4;
                data.specindex2 = He4;
                data.specindex3 = He4;
                break;

            case irg3a:
                data.prefactor  = -1.0_rt;
                data.specindex1 = C12;
                break;

            case iroga:
                data.prefactor  = 1.0_rt;
                data.specindex1 = O16;
                break;

            case ircag:
                data.prefactor  = -1.0_rt;
                data.specindex1 = C12;
                data.specindex2 = He4;
                break;

            case ir1212:
                data.prefactor  = -1.0_rt;
                data.specindex1 = C12;
                data.specindex2 = C12;
                break;

            case ir1216:
                data.prefactor  = -1.0_rt;
                data.specindex1 = C12;
                data.specindex2 = O16;
                break;

            }
            break;

        case O16:

            switch (rateindex) {

            case iroga:
                data.prefactor  = -1.0_rt;
                data.specindex1 = O16;
                break;

            case ircag:
                data.prefactor  = 1.0_rt;
                data.specindex1 = C12;
                data.specindex2 = He4;
                break;

            case ir1216:
                data.prefactor  = -1.0_rt;
                data.specindex1 = C12;
                data.specindex2 = O16;
                break;

            case ir1616:
                data.prefactor  = -1.0_rt;
                data.specindex1 = O16;
                data.specindex2 = O16;
                break;

            case iroag:
                data.prefactor  = -1.0_rt;
                data.specindex1 = O16;
                data.specindex2 = He4;
                break;

            case irnega:
                data.prefactor  = 1.0_rt;
                data.specindex1 = Ne20;
                break;

            }
            break;

        case Ne20:

            switch (rateindex) {

            case ir1212:
                data.prefactor  = 0.5_rt;
                data.specindex1 = C12;
                data.specindex2 = C12;
                break;

            case iroag:
                data.prefactor  = 1.0_rt;
                data.specindex1 = O16;
                data.specindex2 = He4;
                break;

            case irnega:
                data.prefactor  = -1.0_rt;
                data.specindex1 = Ne20;
                break;

            case irmgga:
                data.prefactor  = 1.0_rt;
                data.specindex1 = Mg24;
                break;

            case irneag:
                data.prefactor  = -1.0_rt;
                data.specindex1 = Ne20;
                data.specindex2 = He4;
                break;

            }
            break;

        case Mg24:

            switch (rateindex) {

            case ir1216:
                data.prefactor  = 0.5_rt;
                data.specindex1 = C12;
                data.specindex2 = O16;
                break;

            case irmgga:
                data.prefactor  = -1.0_rt;
                data.specindex1 = Mg24;
                break;

            case irneag:
                data.prefactor  = 1.0_rt;
                data.specindex1 = Ne20;
                data.specindex2 = He4;
                break;

            case irsiga:
                data.prefactor  = 1.0_rt;
                data.specindex1 = Si28;
                break;

            case irmgag:
                data.prefactor  = -1.0_rt;
                data.specindex1 = Mg24;
                data.specindex2 = He4;
                break;

            }
            break;

        case Si28:

            switch (rateindex) {

            case ir1216:
                data.prefactor  = 0.5_rt;
                data.specindex1 = C12;
                data.specindex2 = O16;
                break;

            case ir1616:
                data.prefactor  = 0.5_rt;
                data.specindex1 = O16;
                data.specindex2 = O16;
                break;

            case irsiga:
                data.prefactor  = -1.0_rt;
                data.specindex1 = Si28;
                break;

            case irmgag:
                data.prefactor  = 1.0_rt;
                data.specindex1 = Mg24;
                data.specindex2 = He4;
                break;

            case irsi2ni:
                data.prefactor  = -1.0_rt;
                data.specindex1 = He4;
                break;

            case irni2si:
                data.prefactor  = 1.0_rt;
                data.specindex1 = Ni56;
                break;

            }
            break;

        case Ni56:

            switch (rateindex) {

            case irsi2ni:
                data.prefactor = 1.0_rt;
                data.specindex1 = He4;
                break;

            case irni2si:
                data.prefactor = -1.0_rt;
                data.specindex1 = Ni56;
                break;

            }
            break;

        }

        return data;
    }

}

#endif
