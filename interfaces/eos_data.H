#ifndef _eos_data_H_
#define _eos_data_H_

#include <AMReX.H>
#include <AMReX_REAL.H>

// store some global limits set at initialization time

namespace EOSData
{
  inline bool initialized;
  inline AMREX_GPU_MANAGED amrex::Real mintemp;
  inline AMREX_GPU_MANAGED amrex::Real maxtemp;
  inline AMREX_GPU_MANAGED amrex::Real mindens;
  inline AMREX_GPU_MANAGED amrex::Real maxdens;
  inline AMREX_GPU_MANAGED amrex::Real minx;
  inline AMREX_GPU_MANAGED amrex::Real maxx;
  inline AMREX_GPU_MANAGED amrex::Real minye;
  inline AMREX_GPU_MANAGED amrex::Real maxye;
  inline AMREX_GPU_MANAGED amrex::Real mine;
  inline AMREX_GPU_MANAGED amrex::Real maxe;
  inline AMREX_GPU_MANAGED amrex::Real minp;
  inline AMREX_GPU_MANAGED amrex::Real maxp;
  inline AMREX_GPU_MANAGED amrex::Real mins;
  inline AMREX_GPU_MANAGED amrex::Real maxs;
  inline AMREX_GPU_MANAGED amrex::Real minh;
  inline AMREX_GPU_MANAGED amrex::Real maxh;
}

#endif
