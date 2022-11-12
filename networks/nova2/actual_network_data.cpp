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

    // binding energies per nucleon in MeV
    amrex::Array1D<amrex::Real, 1, NumSpec> ebind_per_nucleon;

    ebind_per_nucleon(H1) = 0.0_rt;
    ebind_per_nucleon(H2) = 1.112283_rt;
    ebind_per_nucleon(He3) = 2.5726799999999996_rt;
    ebind_per_nucleon(He4) = 7.073915_rt;
    ebind_per_nucleon(Be7) = 5.371548_rt;
    ebind_per_nucleon(B8) = 4.717155_rt;
    ebind_per_nucleon(C12) = 7.680144_rt;
    ebind_per_nucleon(C13) = 7.469849_rt;
    ebind_per_nucleon(N13) = 7.238863_rt;
    ebind_per_nucleon(N14) = 7.475613999999999_rt;
    ebind_per_nucleon(N15) = 7.69946_rt;
    ebind_per_nucleon(O14) = 7.052278_rt;
    ebind_per_nucleon(O15) = 7.463692_rt;
    ebind_per_nucleon(O16) = 7.976206_rt;
    ebind_per_nucleon(O17) = 7.7507280000000005_rt;
    ebind_per_nucleon(F17) = 7.542328_rt;
    ebind_per_nucleon(F18) = 7.631638_rt;

    // convert to binding energies per nucleus in MeV
    for (int i = 1; i <= NumSpec; ++i) {
        bion(i) = ebind_per_nucleon(i) * aion[i-1];
    }

    // Set the mass -- this will be in grams
    for (int i = 1; i <= NumSpec; ++i) {
        mion(i) = (aion[i-1] - zion[i-1]) * C::Legacy::m_n + zion[i-1] * (C::Legacy::m_p + C::Legacy::m_e) - bion(i) * C::Legacy::MeV2gr;
    }

}
