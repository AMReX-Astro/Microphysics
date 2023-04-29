#ifndef NSE_INDEX_H
#define NSE_INDEX_H

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <actual_network.H>

namespace NSE_INDEX
{
  extern bool initialized;

  //store photoionization proton
  extern AMREX_GPU_MANAGED int p_index;

  // store n,p,a indices useful for implementing adaptive statistical equilibrium
  extern AMREX_GPU_MANAGED int h1_index;
  extern AMREX_GPU_MANAGED int n_index;
  extern AMREX_GPU_MANAGED int he4_index;

  // array for storing indices to access different rates
  extern AMREX_GPU_MANAGED amrex::Array2D<int, 1, Rates::NumRates, 1, 6> rate_indices;
}

#endif
