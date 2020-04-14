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

namespace Species {
    enum NetworkSpecies {He4 = 1, C12, O16, Ne20, Mg24, Si28, Ni56, NumberSpecies = Ni56};
};

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
};

#endif
