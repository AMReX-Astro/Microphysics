#ifdef REACTIONS
#include <actual_network.H>
#include <actual_rhs.H>
#endif

void network_init()
{
#ifdef REACTIONS
    actual_network_init();
    actual_rhs_init();
#endif
}
