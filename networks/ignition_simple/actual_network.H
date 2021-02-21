#ifndef actual_network_H
#define actual_network_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>

#include <network_properties.H>

using namespace amrex;

namespace C
{
    namespace Legacy
    {
        // These are the values of the constants used in the original ignition_simple
        constexpr amrex::Real clight  = 2.99792458e10_rt;
        constexpr amrex::Real ev2erg  = 1.60217648740e-12_rt;
        constexpr amrex::Real mev2erg = ev2erg*1.0e6_rt;
        constexpr amrex::Real mev2gr  = mev2erg/(clight*clight);

        constexpr amrex::Real mn = 1.67492721184e-24_rt;
        constexpr amrex::Real mp = 1.67262163783e-24_rt;
        constexpr amrex::Real me = 9.1093821545e-28_rt;

        constexpr amrex::Real avo = 6.0221417930e23_rt;
        constexpr amrex::Real enuc_conv2 = -avo*clight*clight;
    }
}

const std::string network_name = "ignition_simple";

namespace ignition_simple
{
    inline AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> bion;
    inline AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> mion;
}

namespace Rates
{
    const int NumRates = 1;
    const int NumGroups = 2;
}

AMREX_INLINE
void actual_network_init ()
{
    using namespace Species;
    using namespace ignition_simple;

    // Binding energies per nucleus in MeV
    bion(C12)  = 92.16294e0_rt;
    bion(O16)  = 127.62093e0_rt;
    bion(Mg24) = 198.2579e0_rt;

    // Set the mass
    for (int i = 1; i <= NumSpec; ++i) {
        mion(i) = (aion[i-1] - zion[i-1]) * C::Legacy::mn + zion[i-1] * (C::Legacy::mp + C::Legacy::me) - bion(i) * C::Legacy::mev2gr;
    }
}

#endif
