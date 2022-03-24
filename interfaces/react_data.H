#ifndef REACT_DATA_H
#define REACT_DATA_H

#include <AMReX.H>
#include <AMReX_REAL.H>

// store some global limits for when we evaluate reactions rates,
// set at initialization time

namespace ReactData
{
    extern bool initialized;
    extern AMREX_GPU_MANAGED amrex::Real mintemp;
    extern AMREX_GPU_MANAGED amrex::Real maxtemp;
}

#endif
