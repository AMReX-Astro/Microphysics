#include <actual_rhs.H>

namespace RateTable
{
    AMREX_GPU_MANAGED RealArray2D<Rates::NRates, nrattab> rattab;
    AMREX_GPU_MANAGED RealArray2D<Rates::NRates, nrattab> drattabdt;
    AMREX_GPU_MANAGED RealArray1D<nrattab> ttab;
}
