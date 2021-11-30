#include <actual_network.H>


namespace network
{
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> bion;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> mion;
}

namespace reaclib_rates
{

    // Temperature coefficient arrays (numbers correspond to reaction
    // numbers in net_info)

    AMREX_GPU_MANAGED amrex::Array2D<amrex::Real, 1, 7, 1, Rates::NumReaclibSets> ctemp_rate;

    // Index into ctemp_rate, dimension 2, where each rate's
    // coefficients start

    AMREX_GPU_MANAGED amrex::Array1D<int, 1, Rates::NrateReaclib> rate_start_idx;

    // Reaction multiplicities-1 (how many rates contribute - 1)

    AMREX_GPU_MANAGED amrex::Array1D<int,  1, Rates::NrateReaclib> rate_extra_mult;

}

void actual_network_init()
{
    using namespace Species;
    using namespace network;

    // binding energies per nucleon in MeV
    amrex::Array1D<amrex::Real, 1, NumSpec> ebind_per_nucleon;

    ebind_per_nucleon(H1) = 0.0_rt;
    ebind_per_nucleon(He4) = 7.073915_rt;
    ebind_per_nucleon(C12) = 7.680144_rt;
    ebind_per_nucleon(C14) = 7.520319000000001_rt;
    ebind_per_nucleon(N13) = 7.238863_rt;
    ebind_per_nucleon(N14) = 7.475613999999999_rt;
    ebind_per_nucleon(O16) = 7.976206_rt;
    ebind_per_nucleon(O18) = 7.767097_rt;
    ebind_per_nucleon(F18) = 7.631638_rt;
    ebind_per_nucleon(Ne20) = 8.03224_rt;
    ebind_per_nucleon(Ne21) = 7.971712999999999_rt;
    ebind_per_nucleon(Mg24) = 8.260709_rt;
    ebind_per_nucleon(Al27) = 8.331553_rt;
    ebind_per_nucleon(Si28) = 8.447744_rt;
    ebind_per_nucleon(P31) = 8.481167_rt;
    ebind_per_nucleon(S32) = 8.493129000000001_rt;
    ebind_per_nucleon(Cl35) = 8.520278000000001_rt;
    ebind_per_nucleon(Ar36) = 8.519909_rt;
    ebind_per_nucleon(K39) = 8.557025_rt;
    ebind_per_nucleon(Ca40) = 8.551303_rt;
    ebind_per_nucleon(Sc43) = 8.530825_rt;
    ebind_per_nucleon(Ti44) = 8.533520000000001_rt;
    ebind_per_nucleon(V47) = 8.582225000000001_rt;
    ebind_per_nucleon(Cr48) = 8.572269_rt;
    ebind_per_nucleon(Mn51) = 8.633772_rt;
    ebind_per_nucleon(Fe52) = 8.609574_rt;
    ebind_per_nucleon(Co55) = 8.669618_rt;
    ebind_per_nucleon(Ni56) = 8.642779_rt;

    // convert to binding energies per nucleus in MeV
    for (int i = 1; i <= NumSpec; ++i) {
        bion(i) = ebind_per_nucleon(i) * aion[i-1];
    }

    // Set the mass -- this will be in grams
    for (int i = 1; i <= NumSpec; ++i) {
        mion(i) = (aion[i-1] - zion[i-1]) * C::Legacy::m_n + zion[i-1] * (C::Legacy::m_p + C::Legacy::m_e) - bion(i) * C::Legacy::MeV2gr;
    }

}
