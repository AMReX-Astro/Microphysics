#ifdef REACTIONS
#include <actual_network.H>
#include <react_data.H>
#ifdef NEW_NETWORK_IMPLEMENTATION
#include <rhs.H>
#else
#include <actual_rhs.H>
#endif
#ifdef NONAKA_PLOT
#include <nonaka_plot.H>
#endif
#endif

#include <extern_parameters.H>

void network_init()
{

#ifdef REACTIONS
    ReactData::mintemp = react_low_cutoff_temp;
    ReactData::maxtemp = MAX_TEMP;
    ReactData::initialized = true;
#endif

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
#endif
}
