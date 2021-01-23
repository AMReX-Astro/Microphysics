#include <AMReX_Vector.H>
#include <actual_network.H>

namespace iso7
{
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> bion;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> mion;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> wion;
}

namespace Rates
{
    amrex::Vector<std::string> names;
}

void actual_network_init()
{
    using namespace Species;
    using namespace iso7;

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

    // set the names of the reaction rates
    {
        using namespace Rates;
        names.resize(NumRates);
        names[C12_He4_to_O16_forward-1]     = "rcag";
        names[C12_He4_to_O16_reverse-1]     = "roga";
        names[He4_He4_He4_to_C12_forward-1] = "r3a";
        names[He4_He4_He4_to_C12_reverse-1] = "rg3a";
        names[C12_C12-1]                    = "r1212";
        names[C12_O16-1]                    = "r1216";
        names[O16_O16-1]                    = "r1616";
        names[O16_He4_to_Ne20_forward-1]    = "roag";
        names[O16_He4_to_Ne20_reverse-1]    = "rnega";
        names[Ne20_He4_to_Mg24_forward-1]   = "rneag";
        names[Ne20_He4_to_Mg24_reverse-1]   = "rmgga";
        names[Mg24_He4_to_Si28_forward-1]   = "rmgag";
        names[Mg24_He4_to_Si28_reverse-1]   = "rsiga";
        names[Ca40_He4_to_Ti44_forward-1]   = "rcaag";
        names[Ca40_He4_to_Ti44_reverse-1]   = "rtiga";
        names[Si28_7He4_to_Ni56_forward-1]  = "rsi2ni";
        names[Si28_7He4_to_Ni56_reverse-1]  = "rni2si";
    }
}
