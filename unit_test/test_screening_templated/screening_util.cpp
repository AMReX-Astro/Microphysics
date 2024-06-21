#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>
#include <AMReX_Loop.H>

#include <variables.H>
#include <network.H>
#include <eos.H>
#include <rhs_type.H>

#include <screen.H>

#include <cmath>

using namespace amrex;
using namespace unit_test_rp;

// for whatever reason, this doesn't work when inlined
template <int do_T_derivatives, typename dual_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void maybe_seed(dual_t& value) {
  if constexpr (do_T_derivatives) {
    autodiff::seed<1>(value, 1.0_rt);
  }
}

void screen_test_C(const Box& bx,
                   const Real dlogrho, const Real dlogT, const Real dmetal,
                   const plot_t& vars,
                   Array4<Real> const sp) {

  const int ih1 = network_spec_index("hydrogen-1");
  if (ih1 < 0) amrex::Error("Error: ih1 not found");

  const int ihe4 = network_spec_index("helium-4");
  if (ihe4 < 0) amrex::Error("Error: ihe4 not found");

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k)
  {

    // set the composition -- approximately solar
    Real metalicity = 0.0 + static_cast<Real> (k) * dmetal;

    Real xn[NumSpec];

    // for now... the screening using 1-based indexing
    Array1D<Real, 1, NumSpec> ymass;

    for (auto& x : xn) {
      x = metalicity / static_cast<Real>(NumSpec - 2);
    }
    xn[ih1] = 0.75_rt - 0.5_rt * metalicity;
    xn[ihe4] = 0.25_rt - 0.5_rt * metalicity;

    for (int n = 0; n < NumSpec; n++) {
      ymass(n+1) = xn[n] / aion[n];
    }

    constexpr int do_T_derivatives = 1;
    using dual_t = std::conditional_t<do_T_derivatives, autodiff::dual, amrex::Real>;
    dual_t temp_zone = std::pow(10.0, std::log10(temp_min) + static_cast<Real>(j)*dlogT);
    maybe_seed<do_T_derivatives>(temp_zone);

    Real dens_zone = std::pow(10.0, std::log10(dens_min) + static_cast<Real>(i)*dlogrho);

    // store default state
    sp(i, j, k, vars.irho) = dens_zone;
    sp(i, j, k, vars.itemp) = static_cast<Real>(temp_zone);
    for (int n = 0; n < NumSpec; n++) {
      sp(i, j, k, vars.ispec+n) = xn[n];
    }

    for (int i = 0; i < unit_test_rp::loops; ++i) {
    plasma_state_t<dual_t> pstate;
    fill_plasma_state(pstate, temp_zone, dens_zone, ymass);

    Real sc1a;
    Real sc1adt = 0;

    constexpr_for<1, Rates::NumRates+1>([&] (auto n) {
      constexpr int rate = n;
      constexpr RHS::rhs_t data = RHS::rhs_data(rate);

      if constexpr (data.screen_forward_reaction == 0 && data.screen_reverse_reaction == 0) {
        return;
      }
      if (vars.iscn(rate).value == -1) {
        return;
      }

      if constexpr (data.exponent_A == 1 && data.exponent_B == 1 && data.exponent_C == 0) {
        // Forward reaction is A + B, screen using these two species

        constexpr amrex::Real Z1 = NetworkProperties::zion(data.species_A);
        constexpr amrex::Real A1 = NetworkProperties::aion(data.species_A);

        constexpr amrex::Real Z2 = NetworkProperties::zion(data.species_B);
        constexpr amrex::Real A2 = NetworkProperties::aion(data.species_B);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        // Require scn_fac to be evaluated at compile time
        static_assert(scn_fac.z1 == Z1);

        actual_screen(pstate, scn_fac, sc1a, sc1adt);
        sp(i, j, k, vars.iscn(rate).value) = sc1a;
        sp(i, j, k, vars.iscn(rate).dt) = sc1adt;
      }

      if constexpr (data.exponent_A == 2 && data.exponent_B == 0 && data.exponent_C == 0) {
        // Forward reaction is A + A, screen using just this species

        constexpr amrex::Real Z1 = NetworkProperties::zion(data.species_A);
        constexpr amrex::Real A1 = NetworkProperties::aion(data.species_A);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z1, A1);

        static_assert(scn_fac.z1 == Z1);

        actual_screen(pstate, scn_fac, sc1a, sc1adt);
        sp(i, j, k, vars.iscn(rate).value) = sc1a;
        sp(i, j, k, vars.iscn(rate).dt) = sc1adt;
      }

      if constexpr (data.exponent_A == 3 && data.exponent_B == 0 && data.exponent_C == 0) {
        // Forward reaction is triple alpha or an equivalent, screen using A + A
        // and then A + X where X has twice the number of protons and neutrons.

        constexpr amrex::Real Z1 = NetworkProperties::zion(data.species_A);
        constexpr amrex::Real A1 = NetworkProperties::aion(data.species_A);

        constexpr auto scn_fac1 = scrn::calculate_screen_factor(Z1, A1, Z1, A1);

        static_assert(scn_fac1.z1 == Z1);

        actual_screen(pstate, scn_fac1, sc1a, sc1adt);
        sp(i, j, k, vars.iscn(rate).value) = sc1a;
        sp(i, j, k, vars.iscn(rate).dt) = sc1adt;

        constexpr amrex::Real Z2 = 2.0_rt * Z1;
        constexpr amrex::Real A2 = 2.0_rt * A1;

        constexpr auto scn_fac2 = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        static_assert(scn_fac2.z1 == Z1);

        actual_screen(pstate, scn_fac2, sc1a, sc1adt);
        sp(i, j, k, vars.iscn(rate).aux_value) = sc1a;
        sp(i, j, k, vars.iscn(rate).aux_dt) = sc1adt;
      }
    });
    }

  });

}
