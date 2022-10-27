#ifndef actual_network_H
#define actual_network_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>

#include <fundamental_constants.H>
#include <network_properties.H>

using namespace amrex;

void actual_network_init();

const std::string network_name = "pynucastro-cxx";

namespace network
{
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> bion;
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> mion;
}

namespace Rates
{

    enum NetworkRates
    {
        k_he4_fe52__ni56 = 1,
        k_p_co55__ni56 = 2,
        k_he4_fe52__p_co55 = 3,
        k_ni56__he4_fe52__derived = 4,
        k_ni56__p_co55__derived = 5,
        k_p_co55__he4_fe52__derived = 6,
        NumRates = k_p_co55__he4_fe52__derived
    };

    // number of reaclib rates

    const int NrateReaclib = 6;

    // number of tabular rates

    const int NrateTabular = 0;

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "he4_fe52__ni56",  // 1,
        "p_co55__ni56",  // 2,
        "he4_fe52__p_co55",  // 3,
        "ni56__he4_fe52__derived",  // 4,
        "ni56__p_co55__derived",  // 5,
        "p_co55__he4_fe52__derived"  // 6,
    };

}

#endif
