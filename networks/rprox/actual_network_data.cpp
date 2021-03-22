#include <AMReX_Vector.H>
#include <actual_network.H>

namespace network
{
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> ebin;
}

void actual_network_init()
{
    using namespace Species;
    using namespace network;

    // Our convention is that binding energy is negative.  The
    // following are the binding energies in MeV.
    ebin(C12)  = 92.16279_rt;
    ebin(O14)  = 98.7325_rt;
    ebin(O15)  = 111.9569_rt;
    ebin(O16)  = 127.6207_rt;
    ebin(F17)  = 128.2211_rt;
    ebin(Mg22) = 168.5768_rt;
    ebin(S30)  = 243.6866_rt;
    ebin(Ni56) = 483.995_rt;
    ebin(He4)  = 28.29599_rt;
    ebin(H1)   = 0.0_rt;

    // convert to erg / g by multiplying by N_A / aion and converting to erg
    for (int i = 1; i <= NumSpec; ++i) {
        ebin(i) = -ebin(i) * C::Legacy::n_A * C::Legacy::MeV2erg / aion[i-1];
    }
}
