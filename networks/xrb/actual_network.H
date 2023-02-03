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
        k_ne18__f18__weak__wc12 = 4,
        k_na20__ne20__weak__wc12 = 5,
        k_na21__ne21__weak__wc12 = 6,
        k_na22__ne22__weak__wc12 = 7,
        k_mg21__na21__weak__wc12 = 8,
        k_mg22__na22__weak__wc12 = 9,
        k_na20__he4_o16__weak__wc12 = 10,
        k_mg21__p_ne20__weak__wc12 = 11,
        k_p_c12__n13 = 12,
        k_he4_c12__o16 = 13,
        k_p_c13__n14 = 14,
        k_p_n13__o14 = 15,
        k_p_n14__o15 = 16,
        k_he4_n14__f18 = 17,
        k_p_n15__o16 = 18,
        k_he4_o14__ne18 = 19,
        k_he4_o15__ne19 = 20,
        k_p_o16__f17 = 21,
        k_he4_o16__ne20 = 22,
        k_p_f17__ne18 = 23,
        k_he4_f17__na21 = 24,
        k_p_f18__ne19 = 25,
        k_he4_f18__na22 = 26,
        k_he4_ne18__mg22 = 27,
        k_p_ne19__na20 = 28,
        k_he4_ne19__mg23 = 29,
        k_p_ne20__na21 = 30,
        k_he4_ne20__mg24 = 31,
        k_p_ne21__na22 = 32,
        k_p_na20__mg21 = 33,
        k_p_na21__mg22 = 34,
        k_p_na22__mg23 = 35,
        k_c12_c12__he4_ne20 = 36,
        k_he4_n13__p_o16 = 37,
        k_p_n15__he4_c12 = 38,
        k_he4_o14__p_f17 = 39,
        k_c12_o16__he4_mg24 = 40,
        k_p_f18__he4_o15 = 41,
        k_he4_f18__p_ne21 = 42,
        k_he4_ne18__p_na21 = 43,
        k_he4_ne19__p_na22 = 44,
        k_p_ne20__he4_f17 = 45,
        k_he4_na20__p_mg23 = 46,
        k_p_mg24__he4_na21 = 47,
        k_he4_he4_he4__c12 = 48,
        NumRates = k_he4_he4_he4__c12
    };

    // number of reaclib rates

    const int NrateReaclib = 48;

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
        "ne18__f18__weak__wc12",  // 4,
        "na20__ne20__weak__wc12",  // 5,
        "na21__ne21__weak__wc12",  // 6,
        "na22__ne22__weak__wc12",  // 7,
        "mg21__na21__weak__wc12",  // 8,
        "mg22__na22__weak__wc12",  // 9,
        "na20__he4_o16__weak__wc12",  // 10,
        "mg21__p_ne20__weak__wc12",  // 11,
        "p_c12__n13",  // 12,
        "he4_c12__o16",  // 13,
        "p_c13__n14",  // 14,
        "p_n13__o14",  // 15,
        "p_n14__o15",  // 16,
        "he4_n14__f18",  // 17,
        "p_n15__o16",  // 18,
        "he4_o14__ne18",  // 19,
        "he4_o15__ne19",  // 20,
        "p_o16__f17",  // 21,
        "he4_o16__ne20",  // 22,
        "p_f17__ne18",  // 23,
        "he4_f17__na21",  // 24,
        "p_f18__ne19",  // 25,
        "he4_f18__na22",  // 26,
        "he4_ne18__mg22",  // 27,
        "p_ne19__na20",  // 28,
        "he4_ne19__mg23",  // 29,
        "p_ne20__na21",  // 30,
        "he4_ne20__mg24",  // 31,
        "p_ne21__na22",  // 32,
        "p_na20__mg21",  // 33,
        "p_na21__mg22",  // 34,
        "p_na22__mg23",  // 35,
        "c12_c12__he4_ne20",  // 36,
        "he4_n13__p_o16",  // 37,
        "p_n15__he4_c12",  // 38,
        "he4_o14__p_f17",  // 39,
        "c12_o16__he4_mg24",  // 40,
        "p_f18__he4_o15",  // 41,
        "he4_f18__p_ne21",  // 42,
        "he4_ne18__p_na21",  // 43,
        "he4_ne19__p_na22",  // 44,
        "p_ne20__he4_f17",  // 45,
        "he4_na20__p_mg23",  // 46,
        "p_mg24__he4_na21",  // 47,
        "he4_he4_he4__c12"  // 48,
    };

}

#endif
