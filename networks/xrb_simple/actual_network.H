#ifndef actual_network_H
#define actual_network_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>

#include <fundamental_constants.H>
#include <network_properties.H>
#include <microphysics_math.H>

using namespace amrex;

void actual_network_init();

const std::string network_name = "xrb-simple";

namespace network
{
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> ebin;
}

namespace Rates
{
    enum NetworkRates
    {
        ir3a    = 1,
        irag15  = 2,
        irap14  = 3,
        irap18  = 4,
        irwk14o = 5,
        irwk15o = 6
    };

    const int NumGroups = 1;

    // hard-coded weak rates from Stan; should double check with newer estimates
    constexpr Real wk14o = 9.832e-3_rt;
    constexpr Real wk15o = 5.682e-3_rt;
};

#endif
