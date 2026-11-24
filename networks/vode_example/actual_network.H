#ifndef actual_network_H
#define actual_network_H

#include <string_view>
#include <AMReX_REAL.H>
#include <AMReX_Array.H>

#include <network_properties.H>

AMREX_INLINE
void actual_network_init() {
}

constexpr std::string_view network_name = "vode_example";

#endif
