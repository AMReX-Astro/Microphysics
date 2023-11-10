#ifndef NSE_TABLE_DATA_H
#define NSE_TABLE_DATA_H

#include <AMReX_Array.H>
#include <AMReX_REAL.H>

using namespace amrex;

#include <network_properties.H>
#include <nse_table_size.H>

namespace nse_table
{

  constexpr int npts = nse_table_size::ntemp * nse_table_size::nden * nse_table_size::nye;

  extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, npts> abartab;
  extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, npts> beatab;
  extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, npts> dyedttab;
  extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, npts> dabardttab;
  extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, npts> dbeadttab;
  extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, npts> enutab;

  extern AMREX_GPU_MANAGED amrex::Array2D<amrex::Real, 1, NumSpec, 1, npts> massfractab;
}

#endif
