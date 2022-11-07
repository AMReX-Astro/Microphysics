#ifndef _nse_index_H
#define _nse_index_H

#include <AMReX.H>
#include <AMReX_REAL.H>

namespace NSE_INDEX
{
  extern bool initialized;

  //store photoionization proton
  extern AMREX_GPU_MANAGED int p_index;

  // store n,p,a indices useful for implementing adaptive statistical equilibrium
  extern AMREX_GPU_MANAGED int h1_index;
  extern AMREX_GPU_MANAGED int n_index;
  extern AMREX_GPU_MANAGED int he4_index;
}

#endif
