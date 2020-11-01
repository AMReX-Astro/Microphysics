#ifndef actual_network_H
#define actual_network_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>

#include <network_properties.H>

using namespace amrex;

void actual_network_init();

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
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> bion;
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> mion;
}

namespace Rates
{
    const int NumRates = 1;
    const int NumGroups = 2;
}

#endif
