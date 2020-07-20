#include <actual_network.H>

namespace ignition_simple
{
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> bion;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> mion;
}

void actual_network_init()
{
    using namespace Species;
    using namespace ignition_simple;

    // Binding energies per nucleus in MeV
    bion(C12)  = 92.16294e0_rt;
    bion(O16)  = 127.62093e0_rt;
    bion(Mg24) = 198.2579e0_rt;

    // Set the mass
    for (int i = 1; i <= NumSpec; ++i) {
        mion(i) = (aion[i-1] - zion[i-1]) * C::Legacy::mn + zion[i-1] * (C::Legacy::mp + C::Legacy::me) - bion(i) * C::Legacy::mev2gr;
    }
}
