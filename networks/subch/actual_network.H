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
        k_he4_he4_he4__c12 = 1,
        k_he4_c12__o16 = 2,
        k_he4_n14__f18 = 3,
        k_he4_f18__p_ne21 = 4,
        k_p_c12__n13 = 5,
        k_he4_n13__p_o16 = 6,
        k_he4_o16__ne20 = 7,
        k_he4_c14__o18 = 8,
        NumRates = k_he4_c14__o18
    };

    // number of reaclib rates

    const int NrateReaclib = 8;
    const int NumReaclibSets = 18;

    // number of tabular rates

    const int NrateTabular = 0;

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "he4_he4_he4__c12",  // 1,
        "he4_c12__o16",  // 2,
        "he4_n14__f18",  // 3,
        "he4_f18__p_ne21",  // 4,
        "p_c12__n13",  // 5,
        "he4_n13__p_o16",  // 6,
        "he4_o16__ne20",  // 7,
        "he4_c14__o18"  // 8,
    };

}

#endif
