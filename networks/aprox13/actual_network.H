#ifndef _actual_network_H_
#define _actual_network_H_

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>

#include <fundamental_constants.H>
#include <network_properties.H>

using namespace amrex;

namespace C
{
    namespace Legacy
    {
        // These are the values of the constants used in the original aprox13
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
    inline AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> bion;
    inline AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> mion;
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
                        irmgap,
                        iralpa,
                        iralpg,
                        irsigp,
                        Si28_He4_to_S32_forward,
                        Si28_He4_to_S32_reverse,
                        irsiap,
                        irppa,
                        irppg,
                        irsgp,
                        S32_He4_to_Ar36_forward,
                        S32_He4_to_Ar36_reverse,
                        irsap,
                        irclpa,
                        irclpg,
                        irargp,
                        Ar36_He4_to_Ca40_forward,
                        Ar36_He4_to_Ca40_reverse,
                        irarap,
                        irkpa,
                        irkpg,
                        ircagp,
                        Ca40_He4_to_Ti44_forward,
                        Ca40_He4_to_Ti44_reverse,
                        ircaap,
                        irscpa,
                        irscpg,
                        irtigp,
                        Ti44_He4_to_Cr48_forward,
                        Ti44_He4_to_Cr48_reverse,
                        irtiap,
                        irvpa,
                        irvpg,
                        ircrgp,
                        Cr48_He4_to_Fe52_forward,
                        Cr48_He4_to_Fe52_reverse,
                        ircrap,
                        irmnpa,
                        irmnpg,
                        irfegp,
                        Fe52_He4_to_Ni56_forward,
                        Fe52_He4_to_Ni56_reverse,
                        irfeap,
                        ircopa,
                        ircopg,
                        irnigp,
                        irs1,
                        irt1,
                        iru1,
                        irv1,
                        irw1,
                        irx1,
                        iry1,
                        NumRates=iry1
                      };

    const int NumGroups = 2;
}

AMREX_INLINE
void actual_network_init ()
{
    using namespace Species;
    using namespace aprox13;

    // Set the binding energy of the element
    bion(He4)  =  28.29603e0_rt;
    bion(C12)  =  92.16294e0_rt;
    bion(O16)  = 127.62093e0_rt;
    bion(Ne20) = 160.64788e0_rt;
    bion(Mg24) = 198.25790e0_rt;
    bion(Si28) = 236.53790e0_rt;
    bion(S32)  = 271.78250e0_rt;
    bion(Ar36) = 306.72020e0_rt;
    bion(Ca40) = 342.05680e0_rt;
    bion(Ti44) = 375.47720e0_rt;
    bion(Cr48) = 411.46900e0_rt;
    bion(Fe52) = 447.70800e0_rt;
    bion(Ni56) = 484.00300e0_rt;

    // Set the mass
    for (int i = 1; i <= NumSpec; ++i) {
        mion(i) = (aion[i-1] - zion[i-1]) * C::Legacy::m_n + zion[i-1] * (C::Legacy::m_p + C::Legacy::m_e) - bion(i) * C::Legacy::MeV2gr;
    }
}

#endif
