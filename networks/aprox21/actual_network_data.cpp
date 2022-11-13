#include <AMReX_Vector.H>
#include <actual_network.H>

namespace network
{
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> bion;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> mion;
}

void actual_network_init()
{
    using namespace Species;
    using namespace network;

    // Set the binding energy of the element
    bion(H1)   = 0.0e0_rt;
    bion(He3)  = 7.71819e0_rt;
    bion(He4)  = 28.29603e0_rt;
    bion(C12)  = 92.16294e0_rt;
    bion(N14)  = 104.65998e0_rt;
    bion(O16)  = 127.62093e0_rt;
    bion(Ne20) = 160.64788e0_rt;
    bion(Mg24) = 198.25790e0_rt;
    bion(Si28) = 236.53790e0_rt;
    bion(S32)  = 271.78250e0_rt;
    bion(Ar36) = 306.72020e0_rt;
    bion(Ca40) = 342.05680e0_rt;
    bion(Ti44) = 375.47720e0_rt;
    bion(Cr48) = 411.46900e0_rt;
    bion(Cr56) = 488.4970e0_rt;
    bion(Fe52) = 447.70800e0_rt;
    bion(Fe54) = 471.7696e0_rt;
    bion(Fe56) = 492.2450e0_rt;
    bion(Ni56) = 484.00300e0_rt;
    bion(N)    = 0.0e0_rt;
    bion(P)    = 0.0e0_rt;

    // Set the mass
    for (int i = 1; i <= NumSpec; ++i) {
        mion(i) = (aion[i-1] - zion[i-1]) * C::Legacy::m_n + zion[i-1] * (C::Legacy::m_p + C::Legacy::m_e) - bion(i) * C::Legacy::MeV2gr;
    }
}
