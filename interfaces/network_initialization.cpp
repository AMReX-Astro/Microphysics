#ifdef REACTIONS
#ifdef NETWORK_HAS_CXX_IMPLEMENTATION
#include <actual_network.H>
#include <actual_rhs.H>
#endif
#endif

void network_init()
{
#ifdef REACTIONS
#ifdef NETWORK_HAS_CXX_IMPLEMENTATION
    actual_network_init();
    actual_rhs_init();
#endif
#endif
}
