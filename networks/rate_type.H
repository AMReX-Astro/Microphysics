#ifndef _rate_type_H_
#define _rate_type_H_

#include <AMReX_REAL.H>
#include <AMReX_Array.H>

#include <actual_network.H>

struct rate_t
{
    // the temperature at which the rates were evaluated
    amrex::Real T_eval;

    // the rate data
    amrex::Array2D<amrex::Real, 1, num_rate_groups, 0, Rates::NRates-1> rates;
};

#endif
