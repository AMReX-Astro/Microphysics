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
        k_n__p__weak__wc12 = 1,
        k_si29__n_si28 = 2,
        k_si30__n_si29 = 3,
        k_p31__p_si30 = 4,
        k_s32__p_p31 = 5,
        k_s32__he4_si28 = 6,
        k_n_si28__si29 = 7,
        k_he4_si28__s32 = 8,
        k_n_si29__si30 = 9,
        k_p_si30__p31 = 10,
        k_p_p31__s32 = 11,
        k_he4_si28__p_p31 = 12,
        k_he4_si29__n_s32 = 13,
        k_p_p31__he4_si28 = 14,
        k_n_s32__he4_si29 = 15,
        NumRates = k_n_s32__he4_si29
    };

    // number of reaclib rates

    const int NrateReaclib = 15;

    // number of tabular rates

    const int NrateTabular = 0;

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "n__p__weak__wc12",  // 1,
        "si29__n_si28",  // 2,
        "si30__n_si29",  // 3,
        "p31__p_si30",  // 4,
        "s32__p_p31",  // 5,
        "s32__he4_si28",  // 6,
        "n_si28__si29",  // 7,
        "he4_si28__s32",  // 8,
        "n_si29__si30",  // 9,
        "p_si30__p31",  // 10,
        "p_p31__s32",  // 11,
        "he4_si28__p_p31",  // 12,
        "he4_si29__n_s32",  // 13,
        "p_p31__he4_si28",  // 14,
        "n_s32__he4_si29"  // 15,
    };

}

#endif
