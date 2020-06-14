#include <actual_rhs.H>

void actual_rhs_init()
{
    rates_init();

    screening_init();

    set_up_screening_factors();
}
