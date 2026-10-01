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
template <int do_T_derivatives, typename number_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void maybe_seed(number_t& value) {
  if constexpr (do_T_derivatives) {
    autodiff::seed(value);
  }
}

// Check independent moment perturbations, rather than just one composition
// direction.  With H and He, dY_H = 2*dM1-dM2 and
// dY_He = (dM2-dM1)/2.  Do not renormalize: these are partial derivatives.
void test_screening_derivatives()
{
    const int ih = network_spec_index("hydrogen-1") + 1;
    const int ihe = network_spec_index("helium-4") + 1;
    AMREX_ALWAYS_ASSERT(ih > 0 && ihe > 0);
    amrex::ParallelFor(12, [=] AMREX_GPU_DEVICE (int k) {
#if SCREEN_METHOD == SCREEN_METHOD_chugunov2009 || SCREEN_METHOD == SCREEN_METHOD_chabrier1998
        // In the weak-coupling limit, cancellation in the free-energy fit
        // amplifies the existing fast_atan derivative approximation.  Use
        // stronger coupling for finite-difference checks of these methods.
        constexpr Real log_rho_min = 8.0_rt;
#else
        constexpr Real log_rho_min = 2.0_rt;
#endif
        const Real rho = std::pow(10.0_rt, log_rho_min + 2.0_rt * (k % 4));
        const Real temp = std::pow(10.0_rt, 7.0_rt + (k / 4));
        Array1D<Real, 1, NumSpec> y{};
        y(ih) = 0.5_rt;
        y(ihe) = 0.125_rt;
        plasma_state_t<screening_dual_t> state;
        fill_plasma_state(state, temp, rho, y);
        constexpr auto pair = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 6.0_rt, 12.0_rt);
        Real h, ht, h1, h2, f, ft, f1, f2;
        actual_log_screen(state, pair, h, ht, h1, h2);
        actual_screen(state, pair, f, ft, f1, f2);

        // The original temperature-only interface must retain its values.
        autodiff::dual td = temp;
        autodiff::seed(td);
        plasma_state_t<autodiff::dual> legacy;
        fill_plasma_state(legacy, td, rho, y);
        Real old_h, old_ht;
        actual_log_screen(legacy, pair, old_h, old_ht);
        AMREX_ALWAYS_ASSERT(std::abs(h-old_h) <= 1.e-12_rt * (1.0_rt+std::abs(h)));
        AMREX_ALWAYS_ASSERT(std::abs(ht-old_ht) <= 1.e-12_rt * (1.0_rt+std::abs(ht)));

        // fast_atan evaluates an approximation but its autodiff rule uses
        // the exact atan derivative.  Allow that existing discrepancy for
        // the two free-energy fits that use it; other methods are tighter.
#if SCREEN_METHOD == SCREEN_METHOD_chugunov2009 || SCREEN_METHOD == SCREEN_METHOD_chabrier1998
        constexpr Real derivative_tol = 5.e-3_rt;
#else
        constexpr Real derivative_tol = 2.e-5_rt;
#endif
        for (int dir = 0; dir < 3; ++dir) {
            const Real step = 1.e-5_rt * (dir == 0 ? temp : 1.0_rt);
            auto yp = y;
            auto ym = y;
            if (dir > 0) {
                const Real dh = dir == 1 ? 2.0_rt * step : -step;
                const Real dhe = dir == 1 ? -0.5_rt * step : 0.5_rt * step;
                yp(ih) += dh; ym(ih) -= dh;
                yp(ihe) += dhe; ym(ihe) -= dhe;
            }
            plasma_state_t<Real> plus, minus;
            fill_plasma_state(plus, temp + (dir == 0 ? step : 0.0_rt), rho, yp);
            fill_plasma_state(minus, temp - (dir == 0 ? step : 0.0_rt), rho, ym);
            const Real hp = actual_log_screen(plus, pair);
            const Real hm = actual_log_screen(minus, pair);
            const Real fp = actual_screen(plus, pair);
            const Real fm = actual_screen(minus, pair);
            const Real hd = dir == 0 ? ht : (dir == 1 ? h1 : h2);
            const Real fd = dir == 0 ? ft : (dir == 1 ? f1 : f2);
            // Compare changes over the perturbation, avoiding division by a
            // nearly zero derivative (e.g. M1 in Debye-Huckel screening).
            AMREX_ALWAYS_ASSERT(std::abs(hp-hm-2.0_rt*step*hd) <=
                                derivative_tol * std::abs(2.0_rt*step*hd) +
                                2.e-12_rt * (1.0_rt+std::abs(hp)+std::abs(hm)));
            AMREX_ALWAYS_ASSERT(std::abs(fp-fm-2.0_rt*step*fd) <=
                                derivative_tol * std::abs(2.0_rt*step*fd) +
                                2.e-12_rt * (1.0_rt+std::abs(fp)+std::abs(fm)));
        }
    });
    amrex::Gpu::synchronize();
    amrex::Print() << "Screening T, M1, M2 derivative checks passed" << std::endl;
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
    using number_t = std::conditional_t<do_T_derivatives, autodiff::dual, amrex::Real>;
    number_t temp_zone = std::pow(10.0, std::log10(temp_min) + static_cast<Real>(j)*dlogT);
    maybe_seed<do_T_derivatives>(temp_zone);

    Real dens_zone = std::pow(10.0, std::log10(dens_min) + static_cast<Real>(i)*dlogrho);

    // store default state
    sp(i, j, k, vars.irho) = dens_zone;
    sp(i, j, k, vars.itemp) = static_cast<Real>(temp_zone);
    for (int n = 0; n < NumSpec; n++) {
      sp(i, j, k, vars.ispec+n) = xn[n];
    }

    for (int loop = 0; loop < unit_test_rp::loops; ++loop) {
    plasma_state_t<number_t> pstate;
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
