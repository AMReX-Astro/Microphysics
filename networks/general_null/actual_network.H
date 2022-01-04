#ifndef actual_network_H
#define actual_network_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>

#include <network_properties.H>

using namespace amrex;

void actual_network_init();

const std::string network_name = "general_null";

namespace network
{}

namespace Rates
{
    const int NumRates = 0;
    const int NumGroups = 0;
}

#endif
