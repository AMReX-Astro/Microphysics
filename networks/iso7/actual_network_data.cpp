#include <AMReX_Vector.H>
#include <actual_network.H>

namespace network
{
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> bion;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> mion;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> wion;
}

void actual_network_init()
{
    using namespace Species;
    using namespace network;

    // Set the binding energy of the element
    bion(He4)  = 28.29603e0_rt;
    bion(C12)  = 92.16294e0_rt;
    bion(O16)  = 127.62093e0_rt;
    bion(Ne20) = 160.64788e0_rt;
    bion(Mg24) = 198.25790e0_rt;
    bion(Si28) = 236.53790e0_rt;
    bion(Ni56) = 484.00300e0_rt;

    // Set the mass
    for (int i = 1; i <= NumSpec; ++i) {
        mion(i) = (aion[i-1] - zion[i-1]) * C::Legacy::m_n + zion[i-1] * (C::Legacy::m_p + C::Legacy::m_e) - bion(i) * C::Legacy::MeV2gr;

        // Molar mass
        wion(i) = C::Legacy::n_A * mion(i);

        // Common approximation to molar mass
        wion(i) = aion[i-1];
    }
}
