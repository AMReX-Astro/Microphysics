#include <actual_network.H>


namespace network
{
}


void actual_network_init()
{

}

void balance_charge(burn_t& state)
{

    // update the number density of electrons due to charge conservation
    state.xn[0] = -state.xn[3] - state.xn[7] + state.xn[1] + state.xn[12] +
                  state.xn[6] + state.xn[4] + state.xn[9] + 2.0 * state.xn[11];

}
