#include <AMReX_Vector.H>
#include <actual_network.H>

namespace aprox21
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
    using namespace aprox21;

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

    // set the names of the reaction rates
    {
        using namespace Rates;
        names.resize(NumRates);

        names[ir3a-1]   = "r_he4_he4_he4_to_c12";
        names[irg3a-1]  = "r_c12_to_he4_he4_he4";
        names[ircag-1]  = "r_c12_ag_o16";
        names[ir1212-1] = "r1212";
        names[ir1216-1] = "r1216";
        names[ir1616-1] = "r1616";
        names[iroga-1]  = "r_o16_ga_c12";
        names[iroag-1]  = "r_o16_ag_ne20";
        names[irnega-1] = "r_ne20_ga_o16";
        names[irneag-1] = "r_ne20_ag_mg24";
        names[irmgga-1] = "r_mg24_ga_ne20";
        names[irmgag-1] = "r_mg24_ag_si28";
        names[irsiga-1] = "r_si28_ga_mg24";
        names[irmgap-1] = "r_mg24_ap_al27";
        names[iralpa-1] = "r_al27_pa_mg24";
        names[iralpg-1] = "r_al27_pg_si28";
        names[irsigp-1] = "r_si28_gp_al27";
        names[irsiag-1] = "r_si28_ag_s32";
        names[irsga-1]  = "r_s32_ga_si28";
        names[irsiap-1] = "r_si28_ap_p31";
        names[irppa-1]  = "r_p31_pa_si28";
        names[irppg-1]  = "r_p31_pg_s32";
        names[irsgp-1]  = "r_s32_gp_p31";
        names[irsag-1]  = "r_s32_ag_ar36";
        names[irarga-1] = "r_ar36_ga_s32";
        names[irsap-1]  = "r_s32_ap_cl35";
        names[irclpa-1] = "r_cl35_pa_s32";
        names[irclpg-1] = "r_cl35_pg_ar36";
        names[irargp-1] = "r_ar36_gp_cl35";
        names[irarag-1] = "r_ar36_ag_ca40";
        names[ircaga-1] = "r_ca40_ga_ar36";
        names[irarap-1] = "r_ar36_ap_k39";
        names[irkpa-1]  = "r_k39_pa_ar36";
        names[irkpg-1]  = "r_k39_pg_ca40";
        names[ircagp-1] = "r_ca40_gp_k39";
        names[ircaag-1] = "r_ca40_ag_ti44";
        names[irtiga-1] = "r_ti44_ga_ca40";
        names[ircaap-1] = "r_ca40_ap_sc43";
        names[irscpa-1] = "r_sc43_pa_ca40";
        names[irscpg-1] = "r_sc43_pg_ti44";
        names[irtigp-1] = "r_ti44_gp_sc43";
        names[irtiag-1] = "r_ti44_ag_cr48";
        names[ircrga-1] = "r_cr48_ga_ti44";
        names[irtiap-1] = "r_ti44_ap_v47";
        names[irvpa-1]  = "r_v47_pa_ti44";
        names[irvpg-1]  = "r_v47_pg_cr48";
        names[ircrgp-1] = "r_cr48_gp_v47";
        names[ircrag-1] = "r_cr48_ag_fe52";
        names[irfega-1] = "r_fe52_ga_cr48";
        names[ircrap-1] = "r_cr48_ap_mn51";
        names[irmnpa-1] = "r_mn51_pa_cr48";
        names[irmnpg-1] = "r_mn51_pg_fe52";
        names[irfegp-1] = "r_fe52_gp_mn51";
        names[irfeag-1] = "r_fe52_ag_ni56";
        names[irniga-1] = "r_ni56_ga_fe52";
        names[irfeap-1] = "r_fe52_ap_co55";
        names[ircopa-1] = "r_co55_pa_fe52";
        names[ircopg-1] = "r_co55_pg_ni56";
        names[irnigp-1] = "r_ni56_gp_co55";

        // for fe54 photodisintegration
        names[ir52ng-1] = "r_fe52_ng_fe53";
        names[ir53gn-1] = "r_fe53_gn_fe52";
        names[ir53ng-1] = "r_fe53_ng_fe54";
        names[ir54gn-1] = "r_fe54_gn_fe53";
        names[irfepg-1] = "r_fe54_pg_co55";
        names[ircogp-1] = "r_co55_gp_fe54";

        // for he4 photodisintegration
        names[irheng-1]  = "r_he3_ng_he4";
        names[irhegn-1]  = "r_he4_gn_he3";
        names[irhng-1]   = "r_h1_ng_h2";
        names[irdgn-1]   = "r_h2_gn_h1";
        names[irdpg-1]   = "r_h2_pg_he3";
        names[irhegp-1]  = "r_he3_gp_h2";

        // for weak reactions
        names[irpen-1]   = "r_prot_to_neut";
        names[irnep-1]   = "r_neut_to_prot";
        names[irn56ec-1] = "r_ni56ec_to_fe56";

        // ppchain
        names[irpp-1]    = "rpp_to_he3";
        names[ir33-1]    = "r_he3_he3_to_h1_h1_he4";
        names[irhe3ag-1] = "r_he3_ag_be7";

        // cno cycles
        names[ircpg-1]   = "r_c12_pg_n13";
        names[irnpg-1]   = "r_n14_pg_o15";
        names[iropg-1]   = "r_o16_pg_f17";
        names[ifa-1]     = "rfa";
        names[ifg-1]     = "rfg";
        names[irnag-1]   = "r_n14_ag_f18";

        // for reactions to fe56 
        names[ir54ng-1]   = "r_fe54_ng_fe55";
        names[ir55gn-1]   = "r_fe55_gn_fe54";
        names[ir55ng-1]   = "r_fe55_ng_fe56";
        names[ir56gn-1]   = "r_fe56_gn_fe55";
        names[irfe54ap-1] = "r_fe54_ap_co57";
        names[irco57pa-1] = "r_co57_pa_fe54";
        names[irfe56pg-1] = "r_fe56_pg_co57";
        names[irco57gp-1] = "r_co57_gp_fe56";


        // the equilibrium links
        names[irr1-1]   = "r_al27_equil";
        names[irs1-1]   = "r_p31_equil";
        names[irt1-1]   = "r_cl35_equil";
        names[iru1-1]   = "r_k39_equil";
        names[irv1-1]   = "r_sc43_equil";
        names[irw1-1]   = "r_v47_equil";
        names[irx1-1]   = "r_mn51_equil";

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

        names[irfe56_aux1-1] = "rfe56aux1";
        names[irfe56_aux2-1] = "rfe56aux2";
        names[irfe56_aux3-1] = "rfe56aux3";
        names[irfe56_aux4-1] = "rfe56aux4";
    }
    
}
