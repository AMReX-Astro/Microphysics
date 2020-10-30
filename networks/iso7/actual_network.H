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

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    constexpr int num_rhs (int species)
    {
        using namespace Species;

        int num_rhs = 0;

        switch (species) {
        case He4:
            num_rhs = 15;
            break;
        case C12:
            num_rhs = 6;
            break;
        case O16:
            num_rhs = 6;
            break;
        case Ne20:
            num_rhs = 5;
            break;
        case Mg24:
            num_rhs = 5;
            break;
        case Si28:
            num_rhs = 6;
            break;
        case Ni56:
            num_rhs = 2;
            break;
        default:
            num_rhs = 0;
            break;
        }

        return num_rhs;
    }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    constexpr amrex::Real prefactor (int species, int index)
    {
        using namespace Species;

        Real prefactor = 1.0_rt;

        switch (species) {
        case He4:
            switch (index) {
            case 1:
                prefactor = 3.0_rt;
                break;
            case 2:
                prefactor = -0.5_rt;
                break;
            case 3:
                prefactor = 1.0_rt;
                break;
            case 4:
                prefactor = -1.0_rt;
                break;
            case 5:
                prefactor = 0.5_rt;
                break;
            case 6:
                prefactor = 0.5_rt;
                break;
            case 7:
                prefactor = 0.5_rt;
                break;
            case 8:
                prefactor = -1.0_rt;
                break;
            case 9:
                prefactor = 1.0_rt;
                break;
            case 10:
                prefactor = 1.0_rt;
                break;
            case 11:
                prefactor = -1.0_rt;
                break;
            case 12:
                prefactor = 1.0_rt;
                break;
            case 13:
                prefactor = -1.0_rt;
                break;
            case 14:
                prefactor = -7.0_rt;
                break;
            case 15:
                prefactor = 7.0_rt;
                break;
            default:
                prefactor = 1.0_rt;
                break;
            }
            break;
        case C12:
            switch (index) {
            case 1:
                prefactor = 1.0_rt / 6.0_rt;
                break;
            case 2:
                prefactor = -1.0_rt;
                break;
            case 3:
                prefactor = 1.0_rt;
                break;
            case 4:
                prefactor = -1.0_rt;
                break;
            case 5:
                prefactor = -1.0_rt;
                break;
            case 6:
                prefactor = -1.0_rt;
                break;
            default:
                prefactor = 1.0_rt;
                break;
            }
            break;
        case O16:
            switch (index) {
            case 1:
                prefactor = -1.0_rt;
                break;
            case 2:
                prefactor = 1.0_rt;
                break;
            case 3:
                prefactor = -1.0_rt;
                break;
            case 4:
                prefactor = -1.0_rt;
                break;
            case 5:
                prefactor = -1.0_rt;
                break;
            case 6:
                prefactor = 1.0_rt;
                break;
            default:
                prefactor = 1.0_rt;
                break;
            }
            break;
        case Ne20:
            switch (index) {
            case 1:
                prefactor = 0.5_rt;
                break;
            case 2:
                prefactor = 1.0_rt;
                break;
            case 3:
                prefactor = -1.0_rt;
                break;
            case 4:
                prefactor = 1.0_rt;
                break;
            case 5:
                prefactor = -1.0_rt;
                break;
            default:
                prefactor = 1.0_rt;
                break;
            }
            break;
        case Mg24:
            switch (index) {
            case 1:
                prefactor = 0.5_rt;
                break;
            case 2:
                prefactor = -1.0_rt;
                break;
            case 3:
                prefactor = 1.0_rt;
                break;
            case 4:
                prefactor = 1.0_rt;
                break;
            case 5:
                prefactor = -1.0_rt;
                break;
            default:
                prefactor = 1.0_rt;
                break;
            }
            break;
        case Si28:
            switch (index) {
            case 1:
                prefactor = 0.5_rt;
                break;
            case 2:
                prefactor = 0.5_rt;
                break;
            case 3:
                prefactor = -1.0_rt;
                break;
            case 4:
                prefactor = 1.0_rt;
                break;
            case 5:
                prefactor = -1.0_rt;
                break;
            case 6:
                prefactor = 1.0_rt;
                break;
            default:
                prefactor = 1.0_rt;
                break;
            }
            break;
        case Ni56:
            switch (index) {
            case 1:
                prefactor = 1.0_rt;
                break;
            case 2:
                prefactor = -1.0_rt;
                break;
            default:
                prefactor = 1.0_rt;
                break;
            }
            break;
        default:
            prefactor = 1.0_rt;
            break;
        }

        return prefactor;
    }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    constexpr int specindex1 (int species, int index)
    {
        using namespace Species;

        int specindex1 = -1;

        switch (species) {
        case He4:
            switch (index) {
            case 1:
                specindex1 = C12;
                break;
            case 2:
                specindex1 = He4;
                break;
            case 3:
                specindex1 = O16;
                break;
            case 4:
                specindex1 = C12;
                break;
            case 5:
                specindex1 = C12;
                break;
            case 6:
                specindex1 = C12;
                break;
            case 7:
                specindex1 = O16;
                break;
            case 8:
                specindex1 = O16;
                break;
            case 9:
                specindex1 = Ne20;
                break;
            case 10:
                specindex1 = Mg24;
                break;
            case 11:
                specindex1 = Ne20;
                break;
            case 12:
                specindex1 = Si28;
                break;
            case 13:
                specindex1 = Mg24;
                break;
            case 14:
                specindex1 = He4;
                break;
            case 15:
                specindex1 = Ni56;
                break;
            default:
                specindex1 = -1;
                break;
            }
            break;
        case C12:
            switch (index) {
            case 1:
                specindex1 = He4;
                break;
            case 2:
                specindex1 = C12;
                break;
            case 3:
                specindex1 = O16;
                break;
            case 4:
                specindex1 = C12;
                break;
            case 5:
                specindex1 = C12;
                break;
            case 6:
                specindex1 = C12;
                break;
            default:
                specindex1 = -1;
                break;
            }
            break;
        case O16:
            switch (index) {
            case 1:
                specindex1 = O16;
                break;
            case 2:
                specindex1 = C12;
                break;
            case 3:
                specindex1 = C12;
                break;
            case 4:
                specindex1 = O16;
                break;
            case 5:
                specindex1 = O16;
                break;
            case 6:
                specindex1 = Ne20;
                break;
            default:
                specindex1 = -1;
                break;
            }
            break;
        case Ne20:
            switch (index) {
            case 1:
                specindex1 = C12;
                break;
            case 2:
                specindex1 = O16;
                break;
            case 3:
                specindex1 = Ne20;
                break;
            case 4:
                specindex1 = Mg24;
                break;
            case 5:
                specindex1 = Ne20;
                break;
            default:
                specindex1 = -1;
                break;
            }
            break;
        case Mg24:
            switch (index) {
            case 1:
                specindex1 = C12;
                break;
            case 2:
                specindex1 = Mg24;
                break;
            case 3:
                specindex1 = Ne20;
                break;
            case 4:
                specindex1 = Si28;
                break;
            case 5:
                specindex1 = Mg24;
                break;
            default:
                specindex1 = -1;
                break;
            }
            break;
        case Si28:
            switch (index) {
            case 1:
                specindex1 = C12;
                break;
            case 2:
                specindex1 = O16;
                break;
            case 3:
                specindex1 = Si28;
                break;
            case 4:
                specindex1 = Mg24;
                break;
            case 5:
                specindex1 = He4;
                break;
            case 6:
                specindex1 = Ni56;
                break;
            default:
                specindex1 = -1;
                break;
            }
            break;
        case Ni56:
            switch (index) {
            case 1:
                specindex1 = He4;
                break;
            case 2:
                specindex1 = Ni56;
                break;
            default:
                specindex1 = -1;
                break;
            }
            break;
        default:
            specindex1 = -1;
            break;
        }

        return specindex1;
    }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    constexpr int specindex2 (int species, int index)
    {
        using namespace Species;

        int specindex2 = -1;

        switch (species) {
        case He4:
            switch (index) {
            case 1:
                specindex2 = -1;
                break;
            case 2:
                specindex2 = He4;
                break;
            case 3:
                specindex2 = -1;
                break;
            case 4:
                specindex2 = He4;
                break;
            case 5:
                specindex2 = C12;
                break;
            case 6:
                specindex2 = O16;
                break;
            case 7:
                specindex2 = O16;
                break;
            case 8:
                specindex2 = He4;
                break;
            case 9:
                specindex2 = -1;
                break;
            case 10:
                specindex2 = -1;
                break;
            case 11:
                specindex2 = He4;
                break;
            case 12:
                specindex2 = -1;
                break;
            case 13:
                specindex2 = He4;
                break;
            case 14:
                specindex2 = -1;
                break;
            case 15:
                specindex2 = -1;
                break;
            default:
                specindex2 = -1;
                break;
            }
            break;
        case C12:
            switch (index) {
            case 1:
                specindex2 = He4;
                break;
            case 2:
                specindex2 = -1;
                break;
            case 3:
                specindex2 = -1;
                break;
            case 4:
                specindex2 = He4;
                break;
            case 5:
                specindex2 = C12;
                break;
            case 6:
                specindex2 = O16;
                break;
            default:
                specindex2 = -1;
                break;
            }
            break;
        case O16:
            switch (index) {
            case 1:
                specindex2 = -1;
                break;
            case 2:
                specindex2 = He4;
                break;
            case 3:
                specindex2 = O16;
                break;
            case 4:
                specindex2 = O16;
                break;
            case 5:
                specindex2 = He4;
                break;
            case 6:
                specindex2 = -1;
                break;
            default:
                specindex2 = -1;
                break;
            }
            break;
        case Ne20:
            switch (index) {
            case 1:
                specindex2 = C12;
                break;
            case 2:
                specindex2 = He4;
                break;
            case 3:
                specindex2 = -1;
                break;
            case 4:
                specindex2 = -1;
                break;
            case 5:
                specindex2 = He4;
                break;
            default:
                specindex2 = -1;
                break;
            }
            break;
        case Mg24:
            switch (index) {
            case 1:
                specindex2 = O16;
                break;
            case 2:
                specindex2 = -1;
                break;
            case 3:
                specindex2 = He4;
                break;
            case 4:
                specindex2 = -1;
                break;
            case 5:
                specindex2 = He4;
                break;
            default:
                specindex2 = -1;
                break;
            }
            break;
        case Si28:
            switch (index) {
            case 1:
                specindex2 = O16;
                break;
            case 2:
                specindex2 = O16;
                break;
            case 3:
                specindex2 = -1;
                break;
            case 4:
                specindex2 = He4;
                break;
            case 5:
                specindex2 = -1;
                break;
            case 6:
                specindex2 = -1;
                break;
            default:
                specindex2 = -1;
                break;
            }
            break;
        case Ni56:
            switch (index) {
            case 1:
                specindex2 = -1;
                break;
            case 2:
                specindex2 = -1;
                break;
            default:
                specindex2 = -1;
                break;
            }
            break;
        default:
            specindex2 = -1;
            break;
        }

        return specindex2;
    }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    constexpr int specindex3 (int species, int index)
    {
        using namespace Species;

        int specindex3 = -1;

        switch (species) {
        case He4:
            switch (index) {
            case 1:
                specindex3 = -1;
                break;
            case 2:
                specindex3 = He4;
                break;
            case 3:
                specindex3 = -1;
                break;
            case 4:
                specindex3 = -1;
                break;
            case 5:
                specindex3 = -1;
                break;
            case 6:
                specindex3 = -1;
                break;
            case 7:
                specindex3 = -1;
                break;
            case 8:
                specindex3 = -1;
                break;
            case 9:
                specindex3 = -1;
                break;
            case 10:
                specindex3 = -1;
                break;
            case 11:
                specindex3 = -1;
                break;
            case 12:
                specindex3 = -1;
                break;
            case 13:
                specindex3 = -1;
                break;
            case 14:
                specindex3 = -1;
                break;
            case 15:
                specindex3 = -1;
                break;
            default:
                specindex3 = -1;
                break;
            }
            break;
        case C12:
            switch (index) {
            case 1:
                specindex3 = He4;
                break;
            case 2:
                specindex3 = -1;
                break;
            case 3:
                specindex3 = -1;
                break;
            case 4:
                specindex3 = -1;
                break;
            case 5:
                specindex3 = -1;
                break;
            case 6:
                specindex3 = -1;
                break;
            default:
                specindex3 = -1;
                break;
            }
            break;
        case O16:
            switch (index) {
            case 1:
                specindex3 = -1;
                break;
            case 2:
                specindex3 = -1;
                break;
            case 3:
                specindex3 = -1;
                break;
            case 4:
                specindex3 = -1;
                break;
            case 5:
                specindex3 = -1;
                break;
            case 6:
                specindex3 = -1;
                break;
            default:
                specindex3 = -1;
                break;
            }
            break;
        case Ne20:
            switch (index) {
            case 1:
                specindex3 = -1;
                break;
            case 2:
                specindex3 = -1;
                break;
            case 3:
                specindex3 = -1;
                break;
            case 4:
                specindex3 = -1;
                break;
            case 5:
                specindex3 = -1;
                break;
            default:
                specindex3 = -1;
                break;
            }
            break;
        case Mg24:
            switch (index) {
            case 1:
                specindex3 = -1;
                break;
            case 2:
                specindex3 = -1;
                break;
            case 3:
                specindex3 = -1;
                break;
            case 4:
                specindex3 = -1;
                break;
            case 5:
                specindex3 = -1;
                break;
            default:
                specindex3 = -1;
                break;
            }
            break;
        case Si28:
            switch (index) {
            case 1:
                specindex3 = -1;
                break;
            case 2:
                specindex3 = -1;
                break;
            case 3:
                specindex3 = -1;
                break;
            case 4:
                specindex3 = -1;
                break;
            case 5:
                specindex3 = -1;
                break;
            case 6:
                specindex3 = -1;
                break;
            default:
                specindex3 = -1;
                break;
            }
            break;
        case Ni56:
            switch (index) {
            case 1:
                specindex3 = -1;
                break;
            case 2:
                specindex3 = -1;
                break;
            default:
                specindex3 = -1;
                break;
            }
            break;
        default:
            specindex3 = -1;
            break;
        }

        return specindex3;
    }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    constexpr int rateindex (int species, int index)
    {
        using namespace Species;
        using namespace Rates;

        int rateindex = -1;

        switch (species) {
        case He4:
            switch (index) {
            case 1:
                rateindex = irg3a;
                break;
            case 2:
                rateindex = ir3a;
                break;
            case 3:
                rateindex = iroga;
                break;
            case 4:
                rateindex = ircag;
                break;
            case 5:
                rateindex = ir1212;
                break;
            case 6:
                rateindex = ir1216;
                break;
            case 7:
                rateindex = ir1616;
                break;
            case 8:
                rateindex = iroag;
                break;
            case 9:
                rateindex = irnega;
                break;
            case 10:
                rateindex = irmgga;
                break;
            case 11:
                rateindex = irneag;
                break;
            case 12:
                rateindex = irsiga;
                break;
            case 13:
                rateindex = irmgag;
                break;
            case 14:
                rateindex = irsi2ni;
                break;
            case 15:
                rateindex = irni2si;
                break;
            default:
                rateindex = -1;
                break;
            }
            break;
        case C12:
            switch (index) {
            case 1:
                rateindex = ir3a;
                break;
            case 2:
                rateindex = irg3a;
                break;
            case 3:
                rateindex = iroga;
                break;
            case 4:
                rateindex = ircag;
                break;
            case 5:
                rateindex = ir1212;
                break;
            case 6:
                rateindex = ir1216;
                break;
            default:
                rateindex = -1;
                break;
            }
            break;
        case O16:
            switch (index) {
            case 1:
                rateindex = iroga;
                break;
            case 2:
                rateindex = ircag;
                break;
            case 3:
                rateindex = ir1216;
                break;
            case 4:
                rateindex = ir1616;
                break;
            case 5:
                rateindex = iroag;
                break;
            case 6:
                rateindex = irnega;
                break;
            default:
                rateindex = -1;
                break;
            }
            break;
        case Ne20:
            switch (index) {
            case 1:
                rateindex = ir1212;
                break;
            case 2:
                rateindex = iroag;
                break;
            case 3:
                rateindex = irnega;
                break;
            case 4:
                rateindex = irmgga;
                break;
            case 5:
                rateindex = irneag;
                break;
            default:
                rateindex = -1;
                break;
            }
            break;
        case Mg24:
            switch (index) {
            case 1:
                rateindex = ir1216;
                break;
            case 2:
                rateindex = irmgga;
                break;
            case 3:
                rateindex = irneag;
                break;
            case 4:
                rateindex = irsiga;
                break;
            case 5:
                rateindex = irmgag;
                break;
            default:
                rateindex = -1;
                break;
            }
            break;
        case Si28:
            switch (index) {
            case 1:
                rateindex = ir1216;
                break;
            case 2:
                rateindex = ir1616;
                break;
            case 3:
                rateindex = irsiga;
                break;
            case 4:
                rateindex = irmgag;
                break;
            case 5:
                rateindex = irsi2ni;
                break;
            case 6:
                rateindex = irni2si;
                break;
            default:
                rateindex = -1;
                break;
            }
            break;
        case Ni56:
            switch (index) {
            case 1:
                rateindex = irsi2ni;
                break;
            case 2:
                rateindex = irni2si;
                break;
            default:
                rateindex = -1;
                break;
            }
            break;
        default:
            rateindex = -1;
            break;
        }

        return rateindex;
    }

}

#endif
