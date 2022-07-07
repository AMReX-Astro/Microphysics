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

    ebind_per_nucleon(N) = 0.0_rt;
    ebind_per_nucleon(H1) = 0.0_rt;
    ebind_per_nucleon(H2) = 1.112283_rt;
    ebind_per_nucleon(He3) = 2.5726799999999996_rt;
    ebind_per_nucleon(He4) = 7.073915_rt;
    ebind_per_nucleon(Li6) = 5.332331_rt;
    ebind_per_nucleon(Li7) = 5.606439_rt;
    ebind_per_nucleon(Be7) = 5.371548_rt;
    ebind_per_nucleon(Be9) = 6.462668_rt;
    ebind_per_nucleon(B8) = 4.717155_rt;
    ebind_per_nucleon(B10) = 6.475083_rt;
    ebind_per_nucleon(B11) = 6.927732_rt;
    ebind_per_nucleon(C12) = 7.680144_rt;
    ebind_per_nucleon(C13) = 7.469849_rt;
    ebind_per_nucleon(C14) = 7.520319000000001_rt;
    ebind_per_nucleon(N13) = 7.238863_rt;
    ebind_per_nucleon(N14) = 7.475613999999999_rt;
    ebind_per_nucleon(N15) = 7.69946_rt;
    ebind_per_nucleon(O14) = 7.052278_rt;
    ebind_per_nucleon(O15) = 7.463692_rt;
    ebind_per_nucleon(O16) = 7.976206_rt;
    ebind_per_nucleon(O17) = 7.7507280000000005_rt;
    ebind_per_nucleon(O18) = 7.767097_rt;
    ebind_per_nucleon(F17) = 7.542328_rt;
    ebind_per_nucleon(F18) = 7.631638_rt;
    ebind_per_nucleon(F19) = 7.779018_rt;
    ebind_per_nucleon(Ne18) = 7.341257_rt;
    ebind_per_nucleon(Ne19) = 7.567343_rt;
    ebind_per_nucleon(Ne20) = 8.03224_rt;
    ebind_per_nucleon(Ne21) = 7.971712999999999_rt;
    ebind_per_nucleon(Ne22) = 8.080465_rt;
    ebind_per_nucleon(Na21) = 7.765547_rt;
    ebind_per_nucleon(Na22) = 7.915667_rt;
    ebind_per_nucleon(Na23) = 8.111493000000001_rt;
    ebind_per_nucleon(Mg23) = 7.901115_rt;
    ebind_per_nucleon(Mg24) = 8.260709_rt;
    ebind_per_nucleon(Mg25) = 8.223502_rt;
    ebind_per_nucleon(Mg26) = 8.333870000000001_rt;
    ebind_per_nucleon(Al25) = 8.021136_rt;
    ebind_per_nucleon(Al26) = 8.149765_rt;
    ebind_per_nucleon(Al27) = 8.331553_rt;
    ebind_per_nucleon(Si28) = 8.447744_rt;
    ebind_per_nucleon(Si29) = 8.448635_rt;
    ebind_per_nucleon(Si30) = 8.520654_rt;
    ebind_per_nucleon(Si31) = 8.458291_rt;
    ebind_per_nucleon(Si32) = 8.481468000000001_rt;
    ebind_per_nucleon(P29) = 8.251236_rt;
    ebind_per_nucleon(P30) = 8.353506_rt;
    ebind_per_nucleon(P31) = 8.481167_rt;
    ebind_per_nucleon(P32) = 8.464120000000001_rt;
    ebind_per_nucleon(P33) = 8.513806_rt;
    ebind_per_nucleon(S32) = 8.493129000000001_rt;
    ebind_per_nucleon(S33) = 8.49763_rt;
    ebind_per_nucleon(S34) = 8.583497999999999_rt;
    ebind_per_nucleon(S35) = 8.53785_rt;
    ebind_per_nucleon(S36) = 8.575389_rt;
    ebind_per_nucleon(Cl33) = 8.304754999999998_rt;
    ebind_per_nucleon(Cl34) = 8.398969999999998_rt;
    ebind_per_nucleon(Cl35) = 8.520278000000001_rt;
    ebind_per_nucleon(Cl36) = 8.521931_rt;
    ebind_per_nucleon(Cl37) = 8.570281000000001_rt;
    ebind_per_nucleon(Ar36) = 8.519909_rt;
    ebind_per_nucleon(Ar37) = 8.527139_rt;
    ebind_per_nucleon(Ar38) = 8.61428_rt;
    ebind_per_nucleon(Ar39) = 8.562598_rt;
    ebind_per_nucleon(Ar40) = 8.595259_rt;
    ebind_per_nucleon(K37) = 8.339846999999999_rt;
    ebind_per_nucleon(K38) = 8.438058000000002_rt;
    ebind_per_nucleon(K39) = 8.557025_rt;
    ebind_per_nucleon(K40) = 8.53809_rt;
    ebind_per_nucleon(K41) = 8.576072_rt;
    ebind_per_nucleon(Ca40) = 8.551303_rt;
    ebind_per_nucleon(Ca41) = 8.546706_rt;
    ebind_per_nucleon(Ca42) = 8.616563_rt;
    ebind_per_nucleon(Ca43) = 8.600663_rt;
    ebind_per_nucleon(Ca44) = 8.658175_rt;
    ebind_per_nucleon(Ca45) = 8.630545_rt;
    ebind_per_nucleon(Ca46) = 8.668979_rt;
    ebind_per_nucleon(Ca47) = 8.639349_rt;
    ebind_per_nucleon(Ca48) = 8.666686_rt;
    ebind_per_nucleon(Sc43) = 8.530825_rt;
    ebind_per_nucleon(Sc44) = 8.557379000000001_rt;
    ebind_per_nucleon(Sc45) = 8.618931_rt;
    ebind_per_nucleon(Sc46) = 8.622012_rt;
    ebind_per_nucleon(Sc47) = 8.66509_rt;
    ebind_per_nucleon(Sc48) = 8.656203999999999_rt;
    ebind_per_nucleon(Sc49) = 8.686256_rt;
    ebind_per_nucleon(Ti44) = 8.533520000000001_rt;
    ebind_per_nucleon(Ti45) = 8.555722_rt;
    ebind_per_nucleon(Ti46) = 8.656450999999999_rt;
    ebind_per_nucleon(Ti47) = 8.661227_rt;
    ebind_per_nucleon(Ti48) = 8.723006_rt;
    ebind_per_nucleon(Ti49) = 8.711157_rt;
    ebind_per_nucleon(Ti50) = 8.755718_rt;
    ebind_per_nucleon(Ti51) = 8.708988_rt;
    ebind_per_nucleon(V46) = 8.48613_rt;
    ebind_per_nucleon(V47) = 8.582225000000001_rt;
    ebind_per_nucleon(V48) = 8.623061_rt;
    ebind_per_nucleon(V49) = 8.682908_rt;
    ebind_per_nucleon(V50) = 8.695917999999999_rt;
    ebind_per_nucleon(V51) = 8.742099_rt;
    ebind_per_nucleon(V52) = 8.714582_rt;
    ebind_per_nucleon(Cr48) = 8.572269_rt;
    ebind_per_nucleon(Cr49) = 8.613290999999998_rt;
    ebind_per_nucleon(Cr50) = 8.701032_rt;
    ebind_per_nucleon(Cr51) = 8.712005_rt;
    ebind_per_nucleon(Cr52) = 8.775989_rt;
    ebind_per_nucleon(Cr53) = 8.760198_rt;
    ebind_per_nucleon(Cr54) = 8.777955_rt;
    ebind_per_nucleon(Mn50) = 8.532696_rt;
    ebind_per_nucleon(Mn51) = 8.633772_rt;
    ebind_per_nucleon(Mn52) = 8.670328999999999_rt;
    ebind_per_nucleon(Mn53) = 8.734174999999999_rt;
    ebind_per_nucleon(Mn54) = 8.737965_rt;
    ebind_per_nucleon(Mn55) = 8.765022_rt;
    ebind_per_nucleon(Fe52) = 8.609574_rt;
    ebind_per_nucleon(Fe53) = 8.648799_rt;
    ebind_per_nucleon(Fe54) = 8.736381999999999_rt;
    ebind_per_nucleon(Fe55) = 8.746595_rt;
    ebind_per_nucleon(Fe56) = 8.790353999999999_rt;
    ebind_per_nucleon(Fe57) = 8.770279_rt;
    ebind_per_nucleon(Fe58) = 8.79225_rt;
    ebind_per_nucleon(Co53) = 8.477658_rt;
    ebind_per_nucleon(Co54) = 8.569217_rt;
    ebind_per_nucleon(Co55) = 8.669618_rt;
    ebind_per_nucleon(Co56) = 8.694835999999999_rt;
    ebind_per_nucleon(Co57) = 8.741882_rt;
    ebind_per_nucleon(Co58) = 8.738968999999999_rt;
    ebind_per_nucleon(Co59) = 8.768035_rt;
    ebind_per_nucleon(Ni56) = 8.642779_rt;
    ebind_per_nucleon(Ni57) = 8.670933000000002_rt;
    ebind_per_nucleon(Ni58) = 8.732059_rt;
    ebind_per_nucleon(Ni59) = 8.736588_rt;
    ebind_per_nucleon(Ni60) = 8.780774_rt;
    ebind_per_nucleon(Ni61) = 8.765025_rt;
    ebind_per_nucleon(Ni62) = 8.794553_rt;
    ebind_per_nucleon(Ni63) = 8.763493_rt;
    ebind_per_nucleon(Ni64) = 8.777460999999999_rt;
    ebind_per_nucleon(Cu57) = 8.503262000000001_rt;
    ebind_per_nucleon(Cu58) = 8.570967000000001_rt;
    ebind_per_nucleon(Cu59) = 8.642_rt;
    ebind_per_nucleon(Cu60) = 8.665602000000002_rt;
    ebind_per_nucleon(Cu61) = 8.715513999999999_rt;
    ebind_per_nucleon(Cu62) = 8.718081_rt;
    ebind_per_nucleon(Cu63) = 8.752138_rt;
    ebind_per_nucleon(Cu64) = 8.739075000000001_rt;
    ebind_per_nucleon(Cu65) = 8.757095999999999_rt;
    ebind_per_nucleon(Zn59) = 8.473777_rt;
    ebind_per_nucleon(Zn60) = 8.58305_rt;
    ebind_per_nucleon(Zn61) = 8.610308999999999_rt;
    ebind_per_nucleon(Zn62) = 8.679343000000001_rt;
    ebind_per_nucleon(Zn63) = 8.686285_rt;
    ebind_per_nucleon(Zn64) = 8.735905_rt;
    ebind_per_nucleon(Zn65) = 8.724264999999999_rt;
    ebind_per_nucleon(Zn66) = 8.759632_rt;
    ebind_per_nucleon(Ga62) = 8.518642_rt;
    ebind_per_nucleon(Ga63) = 8.583926_rt;
    ebind_per_nucleon(Ga64) = 8.611631_rt;
    ebind_per_nucleon(Ge63) = 8.418716_rt;
    ebind_per_nucleon(Ge64) = 8.528823000000001_rt;

    // convert to binding energies per nucleus in MeV
    for (int i = 1; i <= NumSpec; ++i) {
        bion(i) = ebind_per_nucleon(i) * aion[i-1];
    }

    // Set the mass -- this will be in grams
    for (int i = 1; i <= NumSpec; ++i) {
        mion(i) = (aion[i-1] - zion[i-1]) * C::Legacy::m_n + zion[i-1] * (C::Legacy::m_p + C::Legacy::m_e) - bion(i) * C::Legacy::MeV2gr;
    }

}
