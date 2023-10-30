#include <nse_table_data.H>

#ifdef NSE_TABLE
namespace nse_table
{
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, npts> abartab;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, npts> beatab;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, npts> dyedttab;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, npts> dabardttab;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, npts> enutab;

  AMREX_GPU_MANAGED amrex::Array2D<amrex::Real, 1, NumSpec, 1, npts> massfractab;
}
#endif
