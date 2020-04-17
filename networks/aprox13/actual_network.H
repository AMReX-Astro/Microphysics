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

const std::string network_name = "aprox13";

namespace aprox13
{
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> bion;
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> mion;
}

namespace Rates {
    enum NetworkRates {
                        ir3a = 1,
                        irg3a,
                        ircag,
                        iroga,
                        ir1212,
                        ir1216,
                        ir1616,
                        iroag,
                        irnega,
                        irneag,
                        irmgga,
                        irmgag,
                        irsiga,
                        irmgap,
                        iralpa,
                        iralpg,
                        irsigp,
                        irsiag,
                        irsga,
                        irsiap,
                        irppa,
                        irppg,
                        irsgp,
                        irsag,
                        irarga,
                        irsap,
                        irclpa,
                        irclpg,
                        irargp,
                        irarag,
                        ircaga,
                        irarap,
                        irkpa,
                        irkpg,
                        ircagp,
                        ircaag,
                        irtiga,
                        ircaap,
                        irscpa,
                        irscpg,
                        irtigp,
                        irtiag,
                        ircrga,
                        irtiap,
                        irvpa,
                        irvpg,
                        ircrgp,
                        ircrag,
                        irfega,
                        ircrap,
                        irmnpa,
                        irmnpg,
                        irfegp,
                        irfeag,
                        irniga,
                        irfeap,
                        ircopa,
                        ircopg,
                        irnigp,
                        irr1,
                        irs1,
                        irt1,
                        iru1,
                        irv1,
                        irw1,
                        irx1,
                        iry1,
                        NumRates=iry1,
                        // ir3a_stage2 and irg3a_stage2 are helpers
                        // for indexing the second screening stage
                        // of the triple-alpha reaction.
                        ir3a_stage2,
                        irg3a_stage2
                      };

    const int NumGroups = 2;

    extern amrex::Vector<std::string> names;
};

#endif
