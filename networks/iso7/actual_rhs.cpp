#include <actual_rhs.H>

namespace RateTable
{
    AMREX_GPU_MANAGED Array2D<Real, 1, Rates::NumRates, 1, nrattab> rattab;
    AMREX_GPU_MANAGED Array2D<Real, 1, Rates::NumRates, 1, nrattab> drattabdt;
    AMREX_GPU_MANAGED Array1D<Real, 1, nrattab> ttab;
}
