#include <AMReX_Vector.H>
#include <actual_network.H>
#ifdef NSE_TABLE
#include "nse.H"
#endif

namespace aprox19
{
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> bion;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> mion;
}

#ifdef NSE_TABLE
namespace table
{

  AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, npts> ttlog;
  AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, npts> ddlog;
  AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, npts> yetab;

  AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, npts> abartab;
  AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, npts> ebtab;
  AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, npts> wratetab;

  AMREX_GPU_MANAGED amrex::Array2D<amrex::Real, 1, NumSpec, 1, npts> massfractab;

}
#endif

namespace Rates
{
    amrex::Vector<std::string> names;
}

void actual_network_init()
{
    using namespace Species;
    using namespace aprox19;

#ifdef NSE_TABLE
    init_nse();
#endif

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
    bion(Fe52) = 447.70800e0_rt;
    bion(Fe54) = 471.7696e0_rt;
    bion(Ni56) = 484.00300e0_rt;
    bion(N)    = 0.0e0_rt;
    bion(P)    = 0.0e0_rt;

    // Set the mass
    for (int i = 1; i <= NumSpec; ++i) {
        mion(i) = (aion[i-1] - zion[i-1]) * C::Legacy::m_n + zion[i-1] * (C::Legacy::m_p + C::Legacy::m_e) - bion(i) * C::Legacy::MeV2gr;
    }

    // set the names of the reaction rates
    {
        using namespace Rates;
        names.resize(NumRates);

        names[ir3a-1]   = "r3a  ";
        names[irg3a-1]  = "rg3a ";
        names[ircag-1]  = "rcag ";
        names[ir1212-1] = "r1212";
        names[ir1216-1] = "r1216";
        names[ir1616-1] = "r1616";
        names[iroga-1]  = "roga ";
        names[iroag-1]  = "roag ";
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
        names[irsga-1]  = "rsga ";
        names[irsiap-1] = "rsiap";
        names[irppa-1]  = "rppa ";
        names[irppg-1]  = "rppg ";
        names[irsgp-1]  = "rsgp ";
        names[irsag-1]  = "rsag ";
        names[irarga-1] = "rarga";
        names[irsap-1]  = "rsap ";
        names[irclpa-1] = "rclpa";
        names[irclpg-1] = "rclpg";
        names[irargp-1] = "rargp";
        names[irarag-1] = "rarag";
        names[ircaga-1] = "rcaga";
        names[irarap-1] = "rarap";
        names[irkpa-1]  = "rkpa ";
        names[irkpg-1]  = "rkpg ";
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
        names[irvpa-1]  = "rvpa ";
        names[irvpg-1]  = "rvpg ";
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

        // for fe54 photodisintegration
        names[ir52ng-1] = "r52ng";
        names[ir53gn-1] = "r53gn";
        names[ir53ng-1] = "r53ng";
        names[ir54gn-1] = "r54gn";
        names[irfepg-1] = "rfepg";
        names[ircogp-1] = "rcogp";

        // for he4 photodisintegration
        names[irheng-1]  = "rheng";
        names[irhegn-1]  = "rhegn";
        names[irhng-1]   = "rhng ";
        names[irdgn-1]   = "rdgn ";
        names[irdpg-1]   = "rdpg ";
        names[irhegp-1]  = "rhegp";

        // for weak reactions
        names[irpen-1]   = "rpen ";
        names[irnep-1]   = "rnep ";
        names[irn56ec-1] = "r56ec";

        // ppchain
        names[irpp-1]    = "rpp  ";
        names[ir33-1]    = "r33  ";
        names[irhe3ag-1] = "rhe3ag";

        // cno cycles
        names[ircpg-1]   = "rcpg ";
        names[irnpg-1]   = "rnpg ";
        names[iropg-1]   = "ropg ";
        names[ifa-1]     = "rfa  ";
        names[ifg-1]     = "rfg  ";
        names[irnag-1]   = "rnag ";


        // the dummy links
        names[irr1-1]   = "r1   ";
        names[irs1-1]   = "s1   ";
        names[irt1-1]   = "t1   ";
        names[iru1-1]   = "u1   ";
        names[irv1-1]   = "v1   ";
        names[irw1-1]   = "w1   ";
        names[irx1-1]   = "x1   ";

        names[ir1f54-1] = "r1f54";
        names[ir2f54-1] = "r2f54";
        names[ir3f54-1] = "r3f54";
        names[ir4f54-1] = "r4f54";
        names[ir5f54-1] = "r5f54";
        names[ir6f54-1] = "r6f54";
        names[ir7f54-1] = "r7f54";
        names[ir8f54-1] = "r8f54";

        names[iralf1-1] = "ralf1";
        names[iralf2-1] = "ralf2";
    }
    
}
