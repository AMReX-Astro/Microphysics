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
        k_n13__c13__weak__wc12 = 1,
        k_o14__n14__weak__wc12 = 2,
        k_o15__n15__weak__wc12 = 3,
        k_f17__o17__weak__wc12 = 4,
        k_f18__o18__weak__wc12 = 5,
        k_ne18__f18__weak__wc12 = 6,
        k_ne19__f19__weak__wc12 = 7,
        k_p_c12__n13 = 8,
        k_he4_c12__o16 = 9,
        k_p_c13__n14 = 10,
        k_p_n13__o14 = 11,
        k_p_n14__o15 = 12,
        k_he4_n14__f18 = 13,
        k_p_n15__o16 = 14,
        k_he4_n15__f19 = 15,
        k_he4_o14__ne18 = 16,
        k_he4_o15__ne19 = 17,
        k_p_o16__f17 = 18,
        k_he4_o16__ne20 = 19,
        k_p_o17__f18 = 20,
        k_p_o18__f19 = 21,
        k_p_f17__ne18 = 22,
        k_p_f18__ne19 = 23,
        k_p_f19__ne20 = 24,
        k_he4_ne18__mg22 = 25,
        k_he4_ne20__mg24 = 26,
        k_c12_c12__he4_ne20 = 27,
        k_he4_n13__p_o16 = 28,
        k_p_n15__he4_c12 = 29,
        k_he4_o14__p_f17 = 30,
        k_c12_o16__he4_mg24 = 31,
        k_p_o17__he4_n14 = 32,
        k_p_o18__he4_n15 = 33,
        k_p_f18__he4_o15 = 34,
        k_p_f19__he4_o16 = 35,
        k_p_ne20__he4_f17 = 36,
        k_he4_he4_he4__c12 = 37,
        NumRates = k_he4_he4_he4__c12
    };

    // number of reaclib rates

    const int NrateReaclib = 37;

    // number of tabular rates

    const int NrateTabular = 0;

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "n13__c13__weak__wc12",  // 1,
        "o14__n14__weak__wc12",  // 2,
        "o15__n15__weak__wc12",  // 3,
        "f17__o17__weak__wc12",  // 4,
        "f18__o18__weak__wc12",  // 5,
        "ne18__f18__weak__wc12",  // 6,
        "ne19__f19__weak__wc12",  // 7,
        "p_c12__n13",  // 8,
        "he4_c12__o16",  // 9,
        "p_c13__n14",  // 10,
        "p_n13__o14",  // 11,
        "p_n14__o15",  // 12,
        "he4_n14__f18",  // 13,
        "p_n15__o16",  // 14,
        "he4_n15__f19",  // 15,
        "he4_o14__ne18",  // 16,
        "he4_o15__ne19",  // 17,
        "p_o16__f17",  // 18,
        "he4_o16__ne20",  // 19,
        "p_o17__f18",  // 20,
        "p_o18__f19",  // 21,
        "p_f17__ne18",  // 22,
        "p_f18__ne19",  // 23,
        "p_f19__ne20",  // 24,
        "he4_ne18__mg22",  // 25,
        "he4_ne20__mg24",  // 26,
        "c12_c12__he4_ne20",  // 27,
        "he4_n13__p_o16",  // 28,
        "p_n15__he4_c12",  // 29,
        "he4_o14__p_f17",  // 30,
        "c12_o16__he4_mg24",  // 31,
        "p_o17__he4_n14",  // 32,
        "p_o18__he4_n15",  // 33,
        "p_f18__he4_o15",  // 34,
        "p_f19__he4_o16",  // 35,
        "p_ne20__he4_f17",  // 36,
        "he4_he4_he4__c12"  // 37,
    };

}

#endif
