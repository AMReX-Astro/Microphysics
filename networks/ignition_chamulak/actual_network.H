#ifndef actual_network_H
#define actual_network_H

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>

#include <fundamental_constants.H>
#include <network_properties.H>

using namespace amrex;

void actual_network_init();

const std::string network_name = "ignition_chamulak";

namespace network
{
    // M12_chamulak is the effective number of C12 nuclei destroyed per
    // reaction
    const Real M12_chamulak = 2.93e0_rt;
}

namespace Rates
{
    enum NetworkRates {
        NumRates = 1
    };
};

#endif
