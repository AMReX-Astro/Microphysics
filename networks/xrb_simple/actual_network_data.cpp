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
    ebin(H1) = 0.0_rt;
    ebin(He4) = 28.296006_rt;
    ebin(O14) = 98.733369_rt;
    ebin(O15) = 111.956652_rt;
    ebin(Ne18) = 132.144124_rt;
    ebin(Si25) = 187.005541_rt;
    ebin(Fe56) = 492.2450_rt;  // older value, but not important -- this is inert

    // convert to erg / g by multiplying by N_A / aion and converting to erg
    for (int i = 1; i <= NumSpec; ++i) {
        ebin(i) = -ebin(i) * C::Legacy::n_A * C::Legacy::MeV2erg / aion[i-1];
    }
}
