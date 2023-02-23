#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>

#include <variables.H>
#include <network.H>
#include <eos.H>

#include <screen.H>

#include <cmath>

using namespace amrex;

void screen_test_C(const Box& bx,
                   const Real dlogrho, const Real dlogT, const Real dmetal,
                   const plot_t& vars,
                   Array4<Real> const sp) {
  using namespace Species;

  const int ih1 = network_spec_index("hydrogen-1");
  if (ih1 < 0) amrex::Error("Error: ih1 not found");

  const int ihe3 = network_spec_index("helium-3");
  if (ihe3 < 0) amrex::Error("Error: ihe3 not found");

  const int ihe4 = network_spec_index("helium-4");
  if (ihe4 < 0) amrex::Error("Error: ihe4 not found");

  const int ic12 = network_spec_index("carbon-12");
  if (ic12 < 0) amrex::Error("Error: ic12 not found");

  const int in14 = network_spec_index("nitrogen-14");
  if (in14 < 0) amrex::Error("Error: in14 not found");

  const int io16 = network_spec_index("oxygen-16");
  if (io16 < 0) amrex::Error("Error: io16 not found");

  const int ine20 = network_spec_index("neon-20");
  if (ine20 < 0) amrex::Error("Error: ine20 not found");

  const int img24 = network_spec_index("magnesium-24");
  if (img24 < 0) amrex::Error("Error: img24 not found");

  const int isi28 = network_spec_index("silicon-28");
  if (isi28 < 0) amrex::Error("Error: isi28 not found");

  const int is32 = network_spec_index("sulfur-32");
  if (is32 < 0) amrex::Error("Error: is32 not found");

  const int iar36 = network_spec_index("argon-36");    
  if (iar36 < 0) amrex::Error("Error: iar36 not found");

  const int ica40 = network_spec_index("calcium-40");
  if (ica40 < 0) amrex::Error("Error: ica40 not found");

  const int iti44 = network_spec_index("titanium-44");
  if (iti44 < 0) amrex::Error("Error: iti44 not found");

  const int icr48 = network_spec_index("chromium-48");
  if (icr48 < 0) amrex::Error("Error: icr48 not found");

  const int ife52 = network_spec_index("iron-52");
  if (ife52 < 0) amrex::Error("Error: ife52 not found");

  const int ife54 = network_spec_index("iron-54");
  if (ife54 < 0) amrex::Error("Error: ife54 not found");

  const int ife56 = network_spec_index("iron-56");
  if (ife56 < 0) amrex::Error("Error: ife56 not found");

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k)
  {

    // set the composition -- approximately solar
    Real metalicity = 0.0 + static_cast<Real> (k) * dmetal;

    Real xn[NumSpec];

    // for now... the screening using 1-based indexing
    Array1D<Real, 1, NumSpec> ymass;

    for (int n = 0; n < NumSpec; n++) {
      xn[n] = metalicity / static_cast<Real>(NumSpec - 2);
    }
    xn[ih1] = 0.75_rt - 0.5_rt * metalicity;
    xn[ihe4] = 0.25_rt - 0.5_rt * metalicity;

    for (int n = 0; n < NumSpec; n++) {
      ymass(n+1) = xn[n] / aion[n];
    }

    Real temp_zone = std::pow(10.0, std::log10(temp_min) + static_cast<Real>(j)*dlogT);

    Real dens_zone = std::pow(10.0, std::log10(dens_min) + static_cast<Real>(i)*dlogrho);

    // store default state
    sp(i, j, k, vars.irho) = dens_zone;
    sp(i, j, k, vars.itemp) = temp_zone;
    for (int n = 0; n < NumSpec; n++) {
      sp(i, j, k, vars.ispec+n) = xn[n];
    }

    plasma_state_t pstate;
    fill_plasma_state(pstate, temp_zone, dens_zone, ymass);

    Real sc1a;
    Real sc1adt;
    Real sc1add;

    // 3-alpha
    {
        constexpr Real Z1 = NetworkProperties::zion(He4);
        constexpr Real A1 = NetworkProperties::aion(He4);

        constexpr Real Z2 = NetworkProperties::zion(He4);
        constexpr Real A2 = NetworkProperties::aion(He4);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        // Require scn_fac to be evaluated at compile time
        static_assert(scn_fac.z1 == Z1);

        actual_screen(pstate, scn_fac, sc1a, sc1adt);
        sp(i, j, k, vars.iscn_he4_he4) = sc1a;
        sp(i, j, k, vars.iscn_he4_he4_dt) = sc1adt;
    }

    {
        constexpr Real Z1 = NetworkProperties::zion(He4);
        constexpr Real A1 = NetworkProperties::aion(He4);

        constexpr Real Z2 = 4.0_rt;
        constexpr Real A2 = 8.0_rt;

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        static_assert(scn_fac.z1 == Z1);

        actual_screen(pstate, scn_fac, sc1a, sc1adt);
        sp(i, j, k, vars.iscn_he4_be8) = sc1a;
        sp(i, j, k, vars.iscn_he4_be8_dt) = sc1adt;
    }

    // c12(a,g)o16
    {
        constexpr Real Z1 = NetworkProperties::zion(C12);
        constexpr Real A1 = NetworkProperties::aion(C12);

        constexpr Real Z2 = NetworkProperties::zion(He4);
        constexpr Real A2 = NetworkProperties::aion(He4);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        static_assert(scn_fac.z1 == Z1);

        actual_screen(pstate, scn_fac, sc1a, sc1adt);
        sp(i, j, k, vars.iscn_c12_he4) = sc1a;
        sp(i, j, k, vars.iscn_c12_he4_dt) = sc1adt;
    }

    // c12 + c12
    {
        constexpr Real Z1 = NetworkProperties::zion(C12);
        constexpr Real A1 = NetworkProperties::aion(C12);

        constexpr Real Z2 = NetworkProperties::zion(C12);
        constexpr Real A2 = NetworkProperties::aion(C12);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        static_assert(scn_fac.z1 == Z1);

        actual_screen(pstate, scn_fac, sc1a, sc1adt);
        sp(i, j, k, vars.iscn_c12_c12) = sc1a;
        sp(i, j, k, vars.iscn_c12_c12_dt) = sc1adt;
    }

    // c12 + o16
    {
        constexpr Real Z1 = NetworkProperties::zion(C12);
        constexpr Real A1 = NetworkProperties::aion(C12);

        constexpr Real Z2 = NetworkProperties::zion(O16);
        constexpr Real A2 = NetworkProperties::aion(O16);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        static_assert(scn_fac.z1 == Z1);

        actual_screen(pstate, scn_fac, sc1a, sc1adt);
        sp(i, j, k, vars.iscn_c12_o16) = sc1a;
        sp(i, j, k, vars.iscn_c12_o16_dt) = sc1adt;
    }

    // o16 + o16
    {
        constexpr Real Z1 = NetworkProperties::zion(O16);
        constexpr Real A1 = NetworkProperties::aion(O16);

        constexpr Real Z2 = NetworkProperties::zion(O16);
        constexpr Real A2 = NetworkProperties::aion(O16);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        static_assert(scn_fac.z1 == Z1);

        actual_screen(pstate, scn_fac, sc1a, sc1adt);
        sp(i, j, k, vars.iscn_o16_o16) = sc1a;
        sp(i, j, k, vars.iscn_o16_o16_dt) = sc1adt;
    }

    // o16 + he4
    {
        constexpr Real Z1 = NetworkProperties::zion(O16);
        constexpr Real A1 = NetworkProperties::aion(O16);

        constexpr Real Z2 = NetworkProperties::zion(He4);
        constexpr Real A2 = NetworkProperties::aion(He4);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        static_assert(scn_fac.z1 == Z1);

        actual_screen(pstate, scn_fac, sc1a, sc1adt);
        sp(i, j, k, vars.iscn_o16_he4) = sc1a;
        sp(i, j, k, vars.iscn_o16_he4_dt) = sc1adt;
    }

    // ne20(a,g)mg24
    {
        constexpr Real Z1 = NetworkProperties::zion(Ne20);
        constexpr Real A1 = NetworkProperties::aion(Ne20);

        constexpr Real Z2 = NetworkProperties::zion(He4);
        constexpr Real A2 = NetworkProperties::aion(He4);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        static_assert(scn_fac.z1 == Z1);

        actual_screen(pstate, scn_fac, sc1a, sc1adt);
        sp(i, j, k, vars.iscn_ne20_he4) = sc1a;
        sp(i, j, k, vars.iscn_ne20_he4_dt) = sc1adt;
    }

    // mg24(a,g)si28
    {
        constexpr Real Z1 = NetworkProperties::zion(Mg24);
        constexpr Real A1 = NetworkProperties::aion(Mg24);

        constexpr Real Z2 = NetworkProperties::zion(He4);
        constexpr Real A2 = NetworkProperties::aion(He4);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        static_assert(scn_fac.z1 == Z1);

        actual_screen(pstate, scn_fac, sc1a, sc1adt);
        sp(i, j, k, vars.iscn_mg24_he4) = sc1a;
        sp(i, j, k, vars.iscn_mg24_he4_dt) = sc1adt;
    }

    // al27(p,g)si28
    {
        constexpr Real Z1 = 13.0_rt;
        constexpr Real A1 = 27.0_rt;

        constexpr Real Z2 = NetworkProperties::zion(H1);
        constexpr Real A2 = NetworkProperties::aion(H1);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        static_assert(scn_fac.z1 == Z1);

        actual_screen(pstate, scn_fac, sc1a, sc1adt);
        sp(i, j, k, vars.iscn_al27_p) = sc1a;
        sp(i, j, k, vars.iscn_al27_p_dt) = sc1adt;
    }

    // si28 + he4
    {
        constexpr Real Z1 = NetworkProperties::zion(Si28);
        constexpr Real A1 = NetworkProperties::aion(Si28);

        constexpr Real Z2 = NetworkProperties::zion(He4);
        constexpr Real A2 = NetworkProperties::aion(He4);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        static_assert(scn_fac.z1 == Z1);

        actual_screen(pstate, scn_fac, sc1a, sc1adt);
        sp(i, j, k, vars.iscn_si28_he4) = sc1a;
        sp(i, j, k, vars.iscn_si28_he4_dt) = sc1adt;
    }

    // p31(p,g)s32
    {
        constexpr Real Z1 = 15.0_rt;
        constexpr Real A1 = 31.0_rt;

        constexpr Real Z2 = NetworkProperties::zion(H1);
        constexpr Real A2 = NetworkProperties::aion(H1);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        static_assert(scn_fac.z1 == Z1);

        actual_screen(pstate, scn_fac, sc1a, sc1adt);
        sp(i, j, k, vars.iscn_p31_p) = sc1a;
        sp(i, j, k, vars.iscn_p31_p_dt) = sc1adt;
    }

    // s32 to ar36
    {
        constexpr Real Z1 = NetworkProperties::zion(S32);
        constexpr Real A1 = NetworkProperties::aion(S32);

        constexpr Real Z2 = NetworkProperties::zion(He4);
        constexpr Real A2 = NetworkProperties::aion(He4);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        static_assert(scn_fac.z1 == Z1);

        actual_screen(pstate, scn_fac, sc1a, sc1adt);
        sp(i, j, k, vars.iscn_s32_he4) = sc1a;
        sp(i, j, k, vars.iscn_s32_he4_dt) = sc1adt;
    }

    // cl35(p,g)ar36
    {
        constexpr Real Z1 = 17.0_rt;
        constexpr Real A1 = 35.0_rt;

        constexpr Real Z2 = NetworkProperties::zion(H1);
        constexpr Real A2 = NetworkProperties::aion(H1);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        static_assert(scn_fac.z1 == Z1);

        actual_screen(pstate, scn_fac, sc1a, sc1adt);
        sp(i, j, k, vars.iscn_cl35_p) = sc1a;
        sp(i, j, k, vars.iscn_cl35_p_dt) = sc1adt;
    }

    // ar36 to ca40
    {
        constexpr Real Z1 = NetworkProperties::zion(Ar36);
        constexpr Real A1 = NetworkProperties::aion(Ar36);

        constexpr Real Z2 = NetworkProperties::zion(He4);
        constexpr Real A2 = NetworkProperties::aion(He4);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        static_assert(scn_fac.z1 == Z1);

        actual_screen(pstate, scn_fac, sc1a, sc1adt);
        sp(i, j, k, vars.iscn_ar36_he4) = sc1a;
        sp(i, j, k, vars.iscn_ar36_he4_dt) = sc1adt;
    }

    // k39(p,g)ca40
    {
        constexpr Real Z1 = 19.0_rt;
        constexpr Real A1 = 39.0_rt;

        constexpr Real Z2 = NetworkProperties::zion(H1);
        constexpr Real A2 = NetworkProperties::aion(H1);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        static_assert(scn_fac.z1 == Z1);

        actual_screen(pstate, scn_fac, sc1a, sc1adt);
        sp(i, j, k, vars.iscn_k39_p) = sc1a;
        sp(i, j, k, vars.iscn_k39_p_dt) = sc1adt;
    }

    // ca40 to ti44
    {
        constexpr Real Z1 = NetworkProperties::zion(Ca40);
        constexpr Real A1 = NetworkProperties::aion(Ca40);

        constexpr Real Z2 = NetworkProperties::zion(He4);
        constexpr Real A2 = NetworkProperties::aion(He4);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        static_assert(scn_fac.z1 == Z1);

        actual_screen(pstate, scn_fac, sc1a, sc1adt);
        sp(i, j, k, vars.iscn_ca40_he4) = sc1a;
        sp(i, j, k, vars.iscn_ca40_he4_dt) = sc1adt;
    }

    // sc43(p,g)ti44
    {
        constexpr Real Z1 = 21.0_rt;
        constexpr Real A1 = 43.0_rt;

        constexpr Real Z2 = NetworkProperties::zion(H1);
        constexpr Real A2 = NetworkProperties::aion(H1);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        static_assert(scn_fac.z1 == Z1);

        actual_screen(pstate, scn_fac, sc1a, sc1adt);
        sp(i, j, k, vars.iscn_sc43_p) = sc1a;
        sp(i, j, k, vars.iscn_sc43_p_dt) = sc1adt;
    }

    // ti44 to cr48
    {
        constexpr Real Z1 = NetworkProperties::zion(Ti44);
        constexpr Real A1 = NetworkProperties::aion(Ti44);

        constexpr Real Z2 = NetworkProperties::zion(He4);
        constexpr Real A2 = NetworkProperties::aion(He4);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        static_assert(scn_fac.z1 == Z1);

        actual_screen(pstate, scn_fac, sc1a, sc1adt);
        sp(i, j, k, vars.iscn_ti44_he4) = sc1a;
        sp(i, j, k, vars.iscn_ti44_he4_dt) = sc1adt;
    }

    // v47(p,g)cr48
    {
        constexpr Real Z1 = 23.0_rt;
        constexpr Real A1 = 47.0_rt;

        constexpr Real Z2 = NetworkProperties::zion(H1);
        constexpr Real A2 = NetworkProperties::aion(H1);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        static_assert(scn_fac.z1 == Z1);

        actual_screen(pstate, scn_fac, sc1a, sc1adt);
        sp(i, j, k, vars.iscn_v47_p) = sc1a;
        sp(i, j, k, vars.iscn_v47_p_dt) = sc1adt;
    }

    // cr48 to fe52
    {
        constexpr Real Z1 = NetworkProperties::zion(Cr48);
        constexpr Real A1 = NetworkProperties::aion(Cr48);

        constexpr Real Z2 = NetworkProperties::zion(He4);
        constexpr Real A2 = NetworkProperties::aion(He4);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        static_assert(scn_fac.z1 == Z1);

        actual_screen(pstate, scn_fac, sc1a, sc1adt);
        sp(i, j, k, vars.iscn_cr48_he4) = sc1a;
        sp(i, j, k, vars.iscn_cr48_he4_dt) = sc1adt;
    }

    // mn51(p,g)fe52
    {
        constexpr Real Z1 = 25.0_rt;
        constexpr Real A1 = 51.0_rt;

        constexpr Real Z2 = NetworkProperties::zion(H1);
        constexpr Real A2 = NetworkProperties::aion(H1);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        static_assert(scn_fac.z1 == Z1);

        actual_screen(pstate, scn_fac, sc1a, sc1adt);
        sp(i, j, k, vars.iscn_mn51_p) = sc1a;
        sp(i, j, k, vars.iscn_mn51_p_dt) = sc1adt;
    }

    // fe to ni
    {
        constexpr Real Z1 = NetworkProperties::zion(Fe52);
        constexpr Real A1 = NetworkProperties::aion(Fe52);

        constexpr Real Z2 = NetworkProperties::zion(He4);
        constexpr Real A2 = NetworkProperties::aion(He4);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        static_assert(scn_fac.z1 == Z1);

        actual_screen(pstate, scn_fac, sc1a, sc1adt);
        sp(i, j, k, vars.iscn_fe52_he4) = sc1a;
        sp(i, j, k, vars.iscn_fe52_he4_dt) = sc1adt;
    }

    // co55(p,g)ni56
    {
        constexpr Real Z1 = 27.0_rt;
        constexpr Real A1 = 55.0_rt;

        constexpr Real Z2 = NetworkProperties::zion(H1);
        constexpr Real A2 = NetworkProperties::aion(H1);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        static_assert(scn_fac.z1 == Z1);

        actual_screen(pstate, scn_fac, sc1a, sc1adt);
        sp(i, j, k, vars.iscn_co55_p) = sc1a;
        sp(i, j, k, vars.iscn_co55_p_dt) = sc1adt;
    }

    // fe54(p,g)co55
    {
        constexpr Real Z1 = NetworkProperties::zion(Fe54);
        constexpr Real A1 = NetworkProperties::aion(Fe54);

        constexpr Real Z2 = NetworkProperties::zion(H1);
        constexpr Real A2 = NetworkProperties::aion(H1);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        static_assert(scn_fac.z1 == Z1);

        actual_screen(pstate, scn_fac, sc1a, sc1adt);
        sp(i, j, k, vars.iscn_fe54_p) = sc1a;
        sp(i, j, k, vars.iscn_fe54_p_dt) = sc1adt;
    }

    // fe54(a,p)co57
    {
        constexpr Real Z1 = NetworkProperties::zion(Fe54);
        constexpr Real A1 = NetworkProperties::aion(Fe54);

        constexpr Real Z2 = NetworkProperties::zion(He4);
        constexpr Real A2 = NetworkProperties::aion(He4);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        static_assert(scn_fac.z1 == Z1);

        actual_screen(pstate, scn_fac, sc1a, sc1adt);
        sp(i, j, k, vars.iscn_fe54_he4) = sc1a;
        sp(i, j, k, vars.iscn_fe54_he4_dt) = sc1adt;
    }

    // fe56(p,g)co57
    {
        constexpr Real Z1 = NetworkProperties::zion(Fe56);
        constexpr Real A1 = NetworkProperties::aion(Fe56);

        constexpr Real Z2 = NetworkProperties::zion(H1);
        constexpr Real A2 = NetworkProperties::aion(H1);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        static_assert(scn_fac.z1 == Z1);

        actual_screen(pstate, scn_fac, sc1a, sc1adt);
        sp(i, j, k, vars.iscn_fe56_p) = sc1a;
        sp(i, j, k, vars.iscn_fe56_p_dt) = sc1adt;
    }

    // d + p
    {
        constexpr Real Z1 = 1.0_rt;
        constexpr Real A1 = 2.0_rt;

        constexpr Real Z2 = NetworkProperties::zion(H1);
        constexpr Real A2 = NetworkProperties::aion(H1);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        static_assert(scn_fac.z1 == Z1);

        actual_screen(pstate, scn_fac, sc1a, sc1adt);
        sp(i, j, k, vars.iscn_d_p) = sc1a;
        sp(i, j, k, vars.iscn_d_p_dt) = sc1adt;
    }

    // pp
    {
        constexpr Real Z1 = NetworkProperties::zion(H1);
        constexpr Real A1 = NetworkProperties::aion(H1);

        constexpr Real Z2 = NetworkProperties::zion(H1);
        constexpr Real A2 = NetworkProperties::aion(H1);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        static_assert(scn_fac.z1 == Z1);

        actual_screen(pstate, scn_fac, sc1a, sc1adt);
        sp(i, j, k, vars.iscn_p_p) = sc1a;
        sp(i, j, k, vars.iscn_p_p_dt) = sc1adt;
    }

    // he3 + he3
    {
        constexpr Real Z1 = NetworkProperties::zion(He3);
        constexpr Real A1 = NetworkProperties::aion(He3);

        constexpr Real Z2 = NetworkProperties::zion(He3);
        constexpr Real A2 = NetworkProperties::aion(He3);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        static_assert(scn_fac.z1 == Z1);

        actual_screen(pstate, scn_fac, sc1a, sc1adt);
        sp(i, j, k, vars.iscn_he3_he3) = sc1a;
        sp(i, j, k, vars.iscn_he3_he3_dt) = sc1adt;
    }

    // he3 + he4
    {
        constexpr Real Z1 = NetworkProperties::zion(He3);
        constexpr Real A1 = NetworkProperties::aion(He3);

        constexpr Real Z2 = NetworkProperties::zion(He4);
        constexpr Real A2 = NetworkProperties::aion(He4);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        static_assert(scn_fac.z1 == Z1);

        actual_screen(pstate, scn_fac, sc1a, sc1adt);
        sp(i, j, k, vars.iscn_he3_he4) = sc1a;
        sp(i, j, k, vars.iscn_he3_he4_dt) = sc1adt;
    }

    // c12(p,g)n13
    {
        constexpr Real Z1 = NetworkProperties::zion(C12);
        constexpr Real A1 = NetworkProperties::aion(C12);

        constexpr Real Z2 = NetworkProperties::zion(H1);
        constexpr Real A2 = NetworkProperties::aion(H1);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        static_assert(scn_fac.z1 == Z1);

        actual_screen(pstate, scn_fac, sc1a, sc1adt);
        sp(i, j, k, vars.iscn_c12_p) = sc1a;
        sp(i, j, k, vars.iscn_c12_p_dt) = sc1adt;
    }

    // n14(p,g)o15
    {
        constexpr Real Z1 = NetworkProperties::zion(N14);
        constexpr Real A1 = NetworkProperties::aion(N14);

        constexpr Real Z2 = NetworkProperties::zion(H1);
        constexpr Real A2 = NetworkProperties::aion(H1);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        static_assert(scn_fac.z1 == Z1);

        actual_screen(pstate, scn_fac, sc1a, sc1adt);
        sp(i, j, k, vars.iscn_n14_p) = sc1a;
        sp(i, j, k, vars.iscn_n14_p_dt) = sc1adt;
    }

    // o16(p,g)f17
    {
        constexpr Real Z1 = NetworkProperties::zion(O16);
        constexpr Real A1 = NetworkProperties::aion(O16);

        constexpr Real Z2 = NetworkProperties::zion(H1);
        constexpr Real A2 = NetworkProperties::aion(H1);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        static_assert(scn_fac.z1 == Z1);

        actual_screen(pstate, scn_fac, sc1a, sc1adt);
        sp(i, j, k, vars.iscn_o16_p) = sc1a;
        sp(i, j, k, vars.iscn_o16_p_dt) = sc1adt;
    }

    // n14(a,g)f18
    {
        constexpr Real Z1 = NetworkProperties::zion(N14);
        constexpr Real A1 = NetworkProperties::aion(N14);

        constexpr Real Z2 = NetworkProperties::zion(He4);
        constexpr Real A2 = NetworkProperties::aion(He4);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        static_assert(scn_fac.z1 == Z1);

        actual_screen(pstate, scn_fac, sc1a, sc1adt);
        sp(i, j, k, vars.iscn_n14_he4) = sc1a;
        sp(i, j, k, vars.iscn_n14_he4_dt) = sc1adt;
    }

  });

}
