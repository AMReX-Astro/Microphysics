#include <actual_rhs.H>

#include <AMReX_Print.H>

void actual_rhs_init()
{
    rates_init();

    screening_init();

    if (use_tables)
    {
        amrex::Print() << "\nInitializing iso7 rate table\n";
        RHS::tabulate_rates();
    }
}
