#ifndef actual_network_H
#define actual_network_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>

#include <network_properties.H>
#include <microphysics_math.H>

using namespace amrex;

void actual_network_init();

namespace C
{
    namespace Legacy
    {
        // These are the values of the constants used in the original aprox13

        constexpr amrex::Real MeV2erg = 1.60217646e-6_rt;
        constexpr amrex::Real N_A = 6.0221415e23_rt;
    }
}

const std::string network_name = "rprox";

namespace rprox
{
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> ebin;
}

namespace Rates
{
    enum NetworkRates
    {
        irlambCNO = 1,
        irag15o   = 2,
        irr1      = 3,
        irag16o   = 4,
        irpg16o   = 5,
        irpg17f   = 6,
        irgp17f   = 7,
        irlambda2 = 8,
        irap14o   = 9,
        irs1      = 10,
        irlambda1 = 11,
        ir3a      = 12,
        irpg12c   = 13,
        irwk14o   = 14,
        irwk17f   = 15,
        irwk15o   = 16,
        irLweak   = 17,
        irla2     = 18,
        NumRates  = irla2
    };

    const int NumGroups = 3;

    extern amrex::Vector<std::string> names;

    constexpr int dlambCNOdh1   = 1;
    constexpr int drs1dhe4      = 2;
    constexpr int drr1dh1       = 3;
    constexpr int dlambda1dhe4  = 4;
    constexpr int dlambda2dhe4  = 5;
    constexpr int delta1        = 6;
    constexpr int delta2        = 7;
    constexpr int r56eff        = 8;
    constexpr int dr56effdt     = 9;
};

#endif
