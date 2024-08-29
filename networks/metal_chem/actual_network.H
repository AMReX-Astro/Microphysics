#ifndef actual_network_H
#define actual_network_H

#include <fstream>
#include <sstream>
#include <AMReX_REAL.H>
#include <AMReX_Array.H>

#include <fundamental_constants.H>
#include <network_properties.H>
#include <burn_type.H>

using namespace amrex;

void actual_network_init();

const std::string network_name = "pynucastro-cxx";

namespace network
{
    extern AMREX_GPU_MANAGED Array1D<Real, 1, 10> semenov_x;
    extern AMREX_GPU_MANAGED Array1D<Real, 1, 1000> semenov_y;
    extern AMREX_GPU_MANAGED Array2D<Real, 1, 10, 1, 1000> semenov_z;

}

namespace Rates
{

    enum NetworkRates
    {
    };

}


void balance_charge(burn_t& state);

#endif
