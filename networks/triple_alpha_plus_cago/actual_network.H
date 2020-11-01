#ifndef actual_network_H
#define actual_network_H

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>

#include <fundamental_constants.H>
#include <network_properties.H>

using namespace amrex;

void actual_network_init();

const std::string network_name = "triple_alpha_plus_cago";

namespace triple_alpha_plus_cago
{
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> bion;
}

namespace Rates
{
    enum NetworkRates {
        ir3a   = 1,
        ircago = 2,
        NumRates = ircago
    };

    const int NumGroups = 2;

    extern amrex::Vector<std::string> names;
};

#endif
