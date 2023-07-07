#include <AMReX_Vector.H>
#include <actual_network.H>
#ifdef NSE_TABLE
#include <nse_table.H>
#endif

#ifdef NSE_TABLE
namespace table
{

  AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, npts> abartab;
  AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, npts> ebtab;
  AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, npts> wratetab;

  AMREX_GPU_MANAGED amrex::Array2D<amrex::Real, 1, NumSpec, 1, npts> massfractab;

}
#endif

void actual_network_init()
{
#ifdef NSE_TABLE
    init_nse();
#endif
}
