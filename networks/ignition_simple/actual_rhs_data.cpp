#include <actual_rhs.H>

void actual_rhs_init ()
{
    using namespace Species;

    screening_init();

    int jscr = 0;
    add_screening_factor(jscr++, zion[C12-1], aion[C12-1], zion[C12-1], aion[C12-1]);
}
