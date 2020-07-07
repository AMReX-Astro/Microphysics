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

        names[ir3a-1]   = "r3a";
        names[irg3a-1]  = "rg3a";
        names[ircag-1]  = "rcag";
        names[ir1212-1] = "r1212";
        names[ir1216-1] = "r1216";
        names[ir1616-1] = "r1616";
        names[iroga-1]  = "roga";
        names[iroag-1]  = "roag";
        names[irnega-1] = "rnega";
        names[irneag-1] = "rneag";
        names[irmgga-1] = "rmgga";
        names[irmgag-1] = "rmgag";
        names[irsiga-1] = "rsiga";
        names[irmgap-1] = "rmgap";
        names[iralpa-1] = "ralpa";
        names[iralpg-1] = "ralpg";
        names[irsigp-1] = "rsigp";
        names[irsiag-1] = "rsiag";
        names[irsga-1]  = "rsga";
        names[irsiap-1] = "rsiap";
        names[irppa-1]  = "rppa";
        names[irppg-1]  = "rppg";
        names[irsgp-1]  = "rsgp";
        names[irsag-1]  = "rsag";
        names[irarga-1] = "rarga";
        names[irsap-1]  = "rsap";
        names[irclpa-1] = "rclpa";
        names[irclpg-1] = "rclpg";
        names[irargp-1] = "rargp";
        names[irarag-1] = "rarag";
        names[ircaga-1] = "rcaga";
        names[irarap-1] = "rarap";
        names[irkpa-1]  = "rkpa";
        names[irkpg-1]  = "rkpg";
        names[ircagp-1] = "rcagp";
        names[ircaag-1] = "rcaag";
        names[irtiga-1] = "rtiga";
        names[ircaap-1] = "rcaap";
        names[irscpa-1] = "rscpa";
        names[irscpg-1] = "rscpg";
        names[irtigp-1] = "rtigp";
        names[irtiag-1] = "rtiag";
        names[ircrga-1] = "rcrga";
        names[irtiap-1] = "rtiap";
        names[irvpa-1]  = "rvpa";
        names[irvpg-1]  = "rvpg";
        names[ircrgp-1] = "rcrgp";
        names[ircrag-1] = "rcrag";
        names[irfega-1] = "rfega";
        names[ircrap-1] = "rcrap";
        names[irmnpa-1] = "rmnpa";
        names[irmnpg-1] = "rmnpg";
        names[irfegp-1] = "rfegp";
        names[irfeag-1] = "rfeag";
        names[irniga-1] = "rniga";
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
