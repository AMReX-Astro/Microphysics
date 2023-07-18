#include <actual_network.H>


namespace network
{
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> bion;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> mion;
}

#ifdef NSE_NET
namespace NSE_INDEX
{
    AMREX_GPU_MANAGED amrex::Array2D<int, 1, Rates::NumRates, 1, 7, Order::C> rate_indices {
        -1, -1, 5, -1, 1, 2, -1,
        -1, 1, 2, -1, -1, 5, 1,
        -1, 1, 5, -1, -1, 6, -1,
        -1, 1, 6, -1, -1, 8, -1,
        -1, 0, 7, -1, -1, 8, -1,
        -1, 1, 7, -1, -1, 9, -1,
        -1, 1, 8, -1, -1, 10, -1,
        -1, 0, 9, -1, -1, 10, -1,
        -1, 2, 2, -1, 0, 9, -1,
        -1, 2, 2, -1, 1, 8, -1,
        -1, 1, 6, -1, 0, 7, -1,
        -1, 0, 7, -1, 1, 6, 11,
        -1, 1, 8, -1, 0, 9, -1,
        -1, 0, 9, -1, 1, 8, 13,
        -1, -1, 4, -1, -1, 3, -1,
        -1, -1, 5, -1, -1, 4, -1,
        -1, -1, 3, -1, -1, 4, 15,
        -1, -1, 4, -1, -1, 5, 16
    };
    AMREX_GPU_MANAGED bool initialized = false;
}
#endif

void actual_network_init()
{
    using namespace Species;
    using namespace network;

    // binding energies per nucleon in MeV
    amrex::Array1D<amrex::Real, 1, NumSpec> ebind_per_nucleon;

    ebind_per_nucleon(H1) = 0.0_rt;
    ebind_per_nucleon(He4) = 7.073915_rt;
    ebind_per_nucleon(O16) = 7.976206_rt;
    ebind_per_nucleon(O20) = 7.568569999999999_rt;
    ebind_per_nucleon(F20) = 7.720134_rt;
    ebind_per_nucleon(Ne20) = 8.03224_rt;
    ebind_per_nucleon(Mg24) = 8.260709_rt;
    ebind_per_nucleon(Al27) = 8.331553_rt;
    ebind_per_nucleon(Si28) = 8.447744_rt;
    ebind_per_nucleon(P31) = 8.481167_rt;
    ebind_per_nucleon(S32) = 8.493129000000001_rt;

    // convert to binding energies per nucleus in MeV
    for (int i = 1; i <= NumSpec; ++i) {
        bion(i) = ebind_per_nucleon(i) * aion[i-1];
    }

    // Set the mass -- this will be in grams
    for (int i = 1; i <= NumSpec; ++i) {
        mion(i) = (aion[i-1] - zion[i-1]) * C::Legacy::m_n + zion[i-1] * (C::Legacy::m_p + C::Legacy::m_e) - bion(i) * C::Legacy::MeV2gr;
    }

}
