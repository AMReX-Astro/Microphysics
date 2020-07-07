#include <actual_network.H>
#ifdef REACTIONS
#include <actual_rhs.H>
#endif

void network_init()
{
    actual_network_init();
#ifdef REACTIONS
    actual_rhs_init();
#endif
}
