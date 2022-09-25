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
        k_o16__he4_c12 = 1,
        k_ne20__he4_o16 = 2,
        k_c12__he4_he4_he4 = 3,
        k_he4_c12__o16 = 4,
        k_he4_o16__ne20 = 5,
        k_c12_c12__he4_ne20 = 6,
        k_he4_ne20__c12_c12 = 7,
        k_he4_he4_he4__c12 = 8,
        NumRates = k_he4_he4_he4__c12
    };

    // number of reaclib rates

    const int NrateReaclib = 8;

    // number of tabular rates

    const int NrateTabular = 0;

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "o16__he4_c12",  // 1,
        "ne20__he4_o16",  // 2,
        "c12__he4_he4_he4",  // 3,
        "he4_c12__o16",  // 4,
        "he4_o16__ne20",  // 5,
        "c12_c12__he4_ne20",  // 6,
        "he4_ne20__c12_c12",  // 7,
        "he4_he4_he4__c12"  // 8,
    };

}

#endif
