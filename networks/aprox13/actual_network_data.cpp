#include <AMReX_Vector.H>
#include <actual_network.H>

namespace aprox13
{
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> bion;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> mion;
}

namespace Rates
{
    amrex::Vector<std::string> names;
}

void actual_network_init()
{
    using namespace Species;
    using namespace aprox13;

    // Set the binding energy of the element
    bion(He4)  =  28.29603e0_rt;
    bion(C12)  =  92.16294e0_rt;
    bion(O16)  = 127.62093e0_rt;
    bion(Ne20) = 160.64788e0_rt;
    bion(Mg24) = 198.25790e0_rt;
    bion(Si28) = 236.53790e0_rt;
    bion(S32)  = 271.78250e0_rt;
    bion(Ar36) = 306.72020e0_rt;
    bion(Ca40) = 342.05680e0_rt;
    bion(Ti44) = 375.47720e0_rt;
    bion(Cr48) = 411.46900e0_rt;
    bion(Fe52) = 447.70800e0_rt;
    bion(Ni56) = 484.00300e0_rt;

    // Set the mass
    for (int i = 1; i <= NumSpec; ++i) {
        mion(i) = (aion[i-1] - zion[i-1]) * C::Legacy::m_n + zion[i-1] * (C::Legacy::m_p + C::Legacy::m_e) - bion(i) * C::Legacy::MeV2gr;
    }

    // set the names of the reaction rates
    {
        using namespace Rates;
        names.resize(NumRates);

        names[He4_He4_He4_to_C12_forward-1] = "r3a";
        names[He4_He4_He4_to_C12_reverse-1] = "rg3a";
        names[C12_He4_to_O16_forward-1]     = "rcag";
        names[C12_C12_forward-1]            = "r1212";
        names[C12_O16_forward-1]            = "r1216";
        names[O16_O16_forward-1]            = "r1616";
        names[C12_He4_to_O16_reverse-1]     = "roga";
        names[O16_He4_to_Ne20_forward-1]    = "roag";
        names[O16_He4_to_Ne20_reverse-1]    = "rnega";
        names[Ne20_He4_to_Mg24_forward-1]   = "rneag";
        names[Ne20_He4_to_Mg24_reverse-1]   = "rmgga";
        names[Mg24_He4_to_Si28_forward-1]   = "rmgag";
        names[Mg24_He4_to_Si28_reverse-1]   = "rsiga";
        names[irmgap-1] = "rmgap";
        names[iralpa-1] = "ralpa";
        names[iralpg-1] = "ralpg";
        names[irsigp-1] = "rsigp";
        names[Si28_He4_to_S32_forward-1]    = "rsiag";
        names[Si28_He4_to_S32_reverse-1]    = "rsga";
        names[irsiap-1] = "rsiap";
        names[irppa-1]  = "rppa";
        names[irppg-1]  = "rppg";
        names[irsgp-1]  = "rsgp";
        names[S32_He4_to_Ar36_forward-1]    = "rsag";
        names[S32_He4_to_Ar36_reverse-1]    = "rarga";
        names[irsap-1]  = "rsap";
        names[irclpa-1] = "rclpa";
        names[irclpg-1] = "rclpg";
        names[irargp-1] = "rargp";
        names[Ar36_He4_to_Ca40_forward-1]   = "rarag";
        names[Ar36_He4_to_Ca40_reverse-1]   = "rcaga";
        names[irarap-1] = "rarap";
        names[irkpa-1]  = "rkpa";
        names[irkpg-1]  = "rkpg";
        names[ircagp-1] = "rcagp";
        names[Ca40_He4_to_Ti44_forward-1]   = "rcaag";
        names[Ca40_He4_to_Ti44_reverse-1]   = "rtiga";
        names[ircaap-1] = "rcaap";
        names[irscpa-1] = "rscpa";
        names[irscpg-1] = "rscpg";
        names[irtigp-1] = "rtigp";
        names[Ti44_He4_to_Cr48_forward-1]   = "rtiag";
        names[Ti44_He4_to_Cr48_reverse-1]   = "rcrga";
        names[irtiap-1] = "rtiap";
        names[irvpa-1]  = "rvpa";
        names[irvpg-1]  = "rvpg";
        names[ircrgp-1] = "rcrgp";
        names[Cr48_He4_to_Fe52_forward-1]   = "rcrag";
        names[Cr48_He4_to_Fe52_reverse-1]   = "rfega";
        names[ircrap-1] = "rcrap";
        names[irmnpa-1] = "rmnpa";
        names[irmnpg-1] = "rmnpg";
        names[irfegp-1] = "rfegp";
        names[Fe52_He4_to_Ni56_forward-1]   = "rfeag";
        names[Fe52_He4_to_Ni56_reverse-1]   = "rniga";
        names[irfeap-1] = "rfeap";
        names[ircopa-1] = "rcopa";
        names[ircopg-1] = "rcopg";
        names[irnigp-1] = "rnigp";
        names[irr1-1]   = "r1";
        names[irs1-1]   = "s1";
        names[irt1-1]   = "t1";
        names[iru1-1]   = "u1";
        names[irv1-1]   = "v1";
        names[irw1-1]   = "w1";
        names[irx1-1]   = "x1";
        names[iry1-1]   = "y1";
    }
}
