#ifndef NSE_TABLE_DATA_H
#define NSE_TABLE_DATA_H

#include <AMReX_Array.H>
#include <AMReX_REAL.H>

using namespace amrex;

#include <network_properties.H>

namespace nse_table
{
  constexpr int ntemp = 71;
  constexpr int nden = 31;
  constexpr int nye = 21;

  constexpr int npts = ntemp * nden * nye;

  constexpr amrex::Real dlogT = 0.02_rt;
  constexpr amrex::Real dlogrho = 0.10_rt;
  constexpr amrex::Real dye = 0.005_rt;

  extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, npts> abartab;
  extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, npts> ebtab;
  extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, npts> wratetab;

  extern AMREX_GPU_MANAGED amrex::Array2D<amrex::Real, 1, NumSpec, 1, npts> massfractab;
}

#endif
