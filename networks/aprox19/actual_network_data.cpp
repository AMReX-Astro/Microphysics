#include <AMReX_Vector.H>
#include <actual_network.H>

namespace aprox19
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
    using namespace aprox19;

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
    bion(n)    = 0.0e0_rt;
    bion(p)    = 0.0e0_rt;

    // Set the mass
    for (int i = 1; i <= NumSpec; ++i) {
        mion(i) = (aion[i-1] - zion[i-1]) * C::Legacy::m_n + zion[i-1] * (C::Legacy::m_p + C::Legacy::m_e) - bion(i) * C::Legacy::MeV2gr;
    }

    // set the names of the reaction rates
    {
        using namespace Rates;
        names.resize(NumRates);

        ratenames[ir3a]   = "r3a  ";
        ratenames[irg3a]  = "rg3a ";
        ratenames[ircag]  = "rcag ";
        ratenames[ir1212] = "r1212";
        ratenames[ir1216] = "r1216";
        ratenames[ir1616] = "r1616";
        ratenames[iroga]  = "roga ";
        ratenames[iroag]  = "roag ";
        ratenames[irnega] = "rnega";
        ratenames[irneag] = "rneag";
        ratenames[irmgga] = "rmgga";
        ratenames[irmgag] = "rmgag";
        ratenames[irsiga] = "rsiga";
        ratenames[irmgap] = "rmgap";
        ratenames[iralpa] = "ralpa";
        ratenames[iralpg] = "ralpg";
        ratenames[irsigp] = "rsigp";
        ratenames[irsiag] = "rsiag";
        ratenames[irsga]  = "rsga ";
        ratenames[irsiap] = "rsiap";
        ratenames[irppa]  = "rppa ";
        ratenames[irppg]  = "rppg ";
        ratenames[irsgp]  = "rsgp ";
        ratenames[irsag]  = "rsag ";
        ratenames[irarga] = "rarga";
        ratenames[irsap]  = "rsap ";
        ratenames[irclpa] = "rclpa";
        ratenames[irclpg] = "rclpg";
        ratenames[irargp] = "rargp";
        ratenames[irarag] = "rarag";
        ratenames[ircaga] = "rcaga";
        ratenames[irarap] = "rarap";
        ratenames[irkpa]  = "rkpa ";
        ratenames[irkpg]  = "rkpg ";
        ratenames[ircagp] = "rcagp";
        ratenames[ircaag] = "rcaag";
        ratenames[irtiga] = "rtiga";
        ratenames[ircaap] = "rcaap";
        ratenames[irscpa] = "rscpa";
        ratenames[irscpg] = "rscpg";
        ratenames[irtigp] = "rtigp";
        ratenames[irtiag] = "rtiag";
        ratenames[ircrga] = "rcrga";
        ratenames[irtiap] = "rtiap";
        ratenames[irvpa]  = "rvpa ";
        ratenames[irvpg]  = "rvpg ";
        ratenames[ircrgp] = "rcrgp";
        ratenames[ircrag] = "rcrag";
        ratenames[irfega] = "rfega";
        ratenames[ircrap] = "rcrap";
        ratenames[irmnpa] = "rmnpa";
        ratenames[irmnpg] = "rmnpg";
        ratenames[irfegp] = "rfegp";
        ratenames[irfeag] = "rfeag";
        ratenames[irniga] = "rniga";
        ratenames[irfeap] = "rfeap";
        ratenames[ircopa] = "rcopa";
        ratenames[ircopg] = "rcopg";
        ratenames[irnigp] = "rnigp";

        // for fe54 photodisintegration
        ratenames[ir52ng] = "r52ng";
        ratenames[ir53gn] = "r53gn";
        ratenames[ir53ng] = "r53ng";
        ratenames[ir54gn] = "r54gn";
        ratenames[irfepg] = "rfepg";
        ratenames[ircogp] = "rcogp";

        // for he4 photodisintegration
        ratenames[irheng]  = "rheng";
        ratenames[irhegn]  = "rhegn";
        ratenames[irhng]   = "rhng ";
        ratenames[irdgn]   = "rdgn ";
        ratenames[irdpg]   = "rdpg ";
        ratenames[irhegp]  = "rhegp";

        // for weak reactions
        ratenames[irpen]   = "rpen ";
        ratenames[irnep]   = "rnep ";
        ratenames[irn56ec] = "r56ec";

        // ppchain
        ratenames[irpp]    = "rpp  ";
        ratenames[ir33]    = "r33  ";
        ratenames[irhe3ag] = "rhe3ag";

        // cno cycles
        ratenames[ircpg]   = "rcpg ";
        ratenames[irnpg]   = "rnpg ";
        ratenames[iropg]   = "ropg ";
        ratenames[ifa]     = "rfa  ";
        ratenames[ifg]     = "rfg  ";
        ratenames[irnag]   = "rnag ";


        // the dummy links
        ratenames[irr1]   = "r1   ";
        ratenames[irs1]   = "s1   ";
        ratenames[irt1]   = "t1   ";
        ratenames[iru1]   = "u1   ";
        ratenames[irv1]   = "v1   ";
        ratenames[irw1]   = "w1   ";
        ratenames[irx1]   = "x1   ";

        ratenames[ir1f54] = "r1f54";
        ratenames[ir2f54] = "r2f54";
        ratenames[ir3f54] = "r3f54";
        ratenames[ir4f54] = "r4f54";
        ratenames[ir5f54] = "r5f54";
        ratenames[ir6f54] = "r6f54";
        ratenames[ir7f54] = "r7f54";
        ratenames[ir8f54] = "r8f54";

        ratenames[iralf1] = "ralf1";
        ratenames[iralf2] = "ralf2";
    }
    
}
