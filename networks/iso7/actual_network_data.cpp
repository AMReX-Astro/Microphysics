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
        mion(i) = (aion[i-1] - zion[i-1]) * C::m_n + zion[i-1] * (C::m_p + C::m_e) - bion(i) * C::MeV2gr;

        // Molar mass
        wion(i) = C::n_A * mion(i);

        // Common approximation to molar mass
        wion(i) = aion[i-1];
    }

    // set the names of the reaction rates
    {
        using namespace Rates;
        names.resize(NumRates);
        names[ircag-1]   = "rcag";
        names[iroga-1]   = "roga";
        names[ir3a-1]    = "r3a";
        names[irg3a-1]   = "rg3a";    // inverse rate
        names[ir1212-1]  = "r1212";
        names[ir1216-1]  = "r1216";
        names[ir1616-1]  = "r1616";
        names[iroag-1]   = "roag";
        names[irnega-1]  = "rnega";
        names[irneag-1]  = "rneag";
        names[irmgga-1]  = "rmgga";
        names[irmgag-1]  = "rmgag";
        names[irsiga-1]  = "rsiga";
        names[ircaag-1]  = "rcaag";
        names[irtiga-1]  = "rtiga";
        names[irsi2ni-1] = "rsi2ni";
        names[irni2si-1] = "rni2si";
    }
}
