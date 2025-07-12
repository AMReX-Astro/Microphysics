#ifndef actual_network_H
#define actual_network_H

#include <string_view>
#include <AMReX_REAL.H>
#include <AMReX_Array.H>

#include <network_properties.H>

void actual_network_init();

constexpr std::string_view network_name = "general_null";

namespace network
{}

namespace Rates
{
    const int NumRates = 0;
    const int NumGroups = 0;
}

#endif
