#ifdef REACTIONS
#include <actual_network.H>
#ifdef NEW_NETWORK_IMPLEMENTATION
#include <rhs.H>
#else
#include <actual_rhs.H>
#endif
#ifdef NONAKA_PLOT
#include <nonaka_plot.H>
#endif
#ifdef NSE_NET
#include <nse_solver.H>
#endif
#endif

void network_init()
{

#ifdef REACTIONS
#ifdef NONAKA_PLOT
nonaka_init();
#endif
#ifdef NEW_NETWORK_IMPLEMENTATION
    actual_network_init();
    RHS::rhs_init();
#else
    actual_network_init();
    actual_rhs_init();
#endif

#ifdef NSE_NET
    init_nse_net();
#endif

#endif

}
