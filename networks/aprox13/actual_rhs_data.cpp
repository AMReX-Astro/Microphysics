#include <actual_rhs.H>

namespace RateTable
{
    AMREX_GPU_MANAGED Array2D<Real, 1, Rates::NumRates, 1, nrattab> rattab;
    AMREX_GPU_MANAGED Array2D<Real, 1, Rates::NumRates, 1, nrattab> drattabdt;
    AMREX_GPU_MANAGED Array1D<Real, 1, nrattab> ttab;
}

void actual_rhs_init()
{
    rates_init();

    screening_init();

    set_up_screening_factors();

    if (use_tables)
    {
        amrex::Print() << "\nInitializing aprox13 rate table\n";
        set_aprox13rat();
    }
}
