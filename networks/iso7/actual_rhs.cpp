#include <actual_rhs.H>

namespace RateTable
{
    AMREX_GPU_MANAGED RealArray2D<Rates::NumRates, nrattab> rattab;
    AMREX_GPU_MANAGED RealArray2D<Rates::NumRates, nrattab> drattabdt;
    AMREX_GPU_MANAGED RealArray1D<nrattab> ttab;
}
