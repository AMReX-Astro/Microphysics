// Native bit-exact validation of the translated rates.H + screen.H against
// the stock network (host only; included by loom_validate.cpp under
// -DLOOM_VALIDATE_RATES). Network-neutral template; per-port hooks live in
// validate_port.H (PORT lines there, not here).
//
// Over T = 1e4..1e16 K (log grid), every branch temperature from
// validate_port.H +-64 ulps, and 20,000 random points, for every parameter
// variant: tfactors (Real and dual), evaluate_analytical_rate for every rate
// 1..Rates::NumRates (both number types), and the port's direct rate checks.
// Then fill_plasma_state + actual_screen for every screened pair from
// network_data.H (enable_debye_huckel_skip = 0/1) over a rho x T grid and
// compositions built from NumSpec. A sensitivity check at the end (T + 1 ulp
// on the Loom side must change some rate) is reported, not counted.
//
// Snapshot of the triple_alpha_plus_cago port's validator (2026-10-18), made
// generic.
#ifndef VALIDATE_RATES_H
#define VALIDATE_RATES_H

#include <cmath>
#include <random>
#include <string>
#include <vector>

#include <extern_parameters.H>
#include <tfactors.H>
#include <actual_network.H>
#include <screen.H>
#include <react_util.H>
#include <loom_validate.H>

#include "rates.H"
#include "screen.H"
#include <loom_params.H>
#include "validate_port.H"

namespace loomval_rates {

using mpl::Real;
using mpl::Dual;
using loomval_port::cmp_dual;

#define LV_TF_FIELDS(X) X(temp) X(t9) X(t92) X(t93) X(t95) X(t912) X(t932) X(t952) X(t972) \
    X(t913) X(t923) X(t943) X(t953) X(t9i) X(t9i2) X(t9i12) X(t9i32) X(t9i13) X(t9i23)       \
    X(t9i43) X(t9i53) X(lnt9)

inline void cmp_tf(LoomCheck& c, const tf_t<Real>& s, const mpl::tf_t<Real>& l)
{
#define LV_TF(f) loom_check_bits(c, "tf_r." #f, s.f, l.f);
    LV_TF_FIELDS(LV_TF)
#undef LV_TF
}

inline void cmp_tf_d(LoomCheck& c, const tf_t<autodiff::dual>& s, const mpl::tf_t<Dual>& l)
{
#define LV_TF(f) cmp_dual(c, "tf_d." #f, s.f, l.f);
    LV_TF_FIELDS(LV_TF)
#undef LV_TF
}

// stock RHS::evaluate_analytical_rate<rate> vs mpl::evaluate_analytical_rate_{r,d}<rate>
template <int rate>
inline void cmp_dispatch(LoomCheck& c, Real T, Real T_loom, const mpl::RatesParams* rp)
{
    RHS::rhs_state_t<Real> sr{};
    sr.tf = get_tfactors(T);
    RHS::rate_t<Real> ar{0.0, 0.0};
    RHS::evaluate_analytical_rate<rate>(sr, ar);
    const mpl::tf_t<Real> lt = mpl::get_tfactors_r(T_loom);
    const mpl::rate_t<Real> lr = mpl::evaluate_analytical_rate_r<rate>(&lt, rp);
    loom_check_bits(c, "eval_r.fr", ar.fr, lr.fr);
    loom_check_bits(c, "eval_r.rr", ar.rr, lr.rr);

    autodiff::dual td = T;
    autodiff::seed(td);
    RHS::rhs_state_t<autodiff::dual> sd{};
    sd.tf = get_tfactors(td);
    RHS::rate_t<autodiff::dual> ad;
    ad.fr = 0.0_rt;
    ad.rr = 0.0_rt;
    RHS::evaluate_analytical_rate<rate>(sd, ad);
    const mpl::tf_t<Dual> ltd = mpl::get_tfactors_d(mpl::dual_seed(T_loom));
    const mpl::rate_t<Dual> ld = mpl::evaluate_analytical_rate_d<rate>(&ltd, rp);
    cmp_dual(c, "eval_d.fr", ad.fr, ld.fr);
    cmp_dual(c, "eval_d.rr", ad.rr, ld.rr);
}

template <int rate>
struct DispatchAll {
    static void run(LoomCheck& c, Real T, Real T_loom, const mpl::RatesParams* rp)
    {
        if constexpr (rate <= mpl::Rates::NumRates) {
            cmp_dispatch<rate>(c, T, T_loom, rp);
            DispatchAll<rate + 1>::run(c, T, T_loom, rp);
        }
    }
};

// noinline: keeps g++ from constant-folding libm calls (MPFR, correctly
// rounded) on the stock side only, for loop-peeled constant T
// (references/footguns.md, "GCC constant-folds libm calls").
__attribute__((noinline)) void one_T(LoomCheck& c, Real T)
{
    mpl::RatesParams rp{};
    fill_rates_params(&rp);

    const tf_t<Real> st = get_tfactors(T);
    const mpl::tf_t<Real> lt = mpl::get_tfactors_r(T);
    cmp_tf(c, st, lt);
    autodiff::dual td = T;
    autodiff::seed(td);
    const tf_t<autodiff::dual> std_ = get_tfactors(td);
    const mpl::tf_t<Dual> ltd = mpl::get_tfactors_d(mpl::dual_seed(T));
    cmp_tf_d(c, std_, ltd);

    DispatchAll<1>::run(c, T, T, &rp);
    loomval_port::direct_rate_checks(c, st, std_, &lt, &ltd, &rp);
}

__attribute__((noinline)) void screen_one(LoomCheck& c, Real T, Real rho, const std::vector<Real>& X)
{
    mpl::RatesParams rp{};
    fill_rates_params(&rp);

    amrex::Array1D<Real, 1, NumSpec> ys;
    Real yl[NumSpec + 1] = {};
    for (int n = 1; n <= NumSpec; ++n) {
        ys(n) = X[n - 1] * aion_inv[n - 1];
        yl[n] = ys(n);
    }

    plasma_state_t<Real> ps;
    fill_plasma_state(ps, T, rho, ys);
    mpl::plasma_state_t<Real> pl{};
    mpl::fill_plasma_state_r(&pl, T, rho, yl);
#define LV_PS_FIELDS(X) X(qlam0z) X(taufac) X(aa) X(temp) X(zbar) X(z2bar) X(n_e) X(gamma_e_fac)
#define LV_PS(f) loom_check_bits(c, "pstate_r." #f, ps.f, pl.f);
    LV_PS_FIELDS(LV_PS)
#undef LV_PS

    autodiff::dual td = T;
    autodiff::seed(td);
    plasma_state_t<autodiff::dual> psd;
    fill_plasma_state(psd, td, rho, ys);
    mpl::plasma_state_t<Dual> pld{};
    mpl::fill_plasma_state_d(&pld, mpl::dual_seed(T), rho, yl);
#define LV_PS(f) cmp_dual(c, "pstate_d." #f, psd.f, pld.f);
    LV_PS_FIELDS(LV_PS)
#undef LV_PS
#undef LV_PS_FIELDS

    for (int pair = 0; pair < mpl::net_num_screen_pairs; ++pair) {
        const scrn::screen_factors_t sf = scrn::calculate_screen_factor(
            mpl::net_screen_z1(pair), mpl::net_screen_a1(pair),
            mpl::net_screen_z2(pair), mpl::net_screen_a2(pair));
        const mpl::screen_factors_t lf{mpl::net_screen_z1(pair), mpl::net_screen_z2(pair),
            mpl::net_screen_a1(pair), mpl::net_screen_a2(pair), mpl::net_screen_zs53(pair),
            mpl::net_screen_z1_53(pair), mpl::net_screen_z2_53(pair), mpl::net_screen_zs13(pair),
            mpl::net_screen_zs13inv(pair), mpl::net_screen_zhat(pair), mpl::net_screen_zhat2(pair),
            mpl::net_screen_lzav(pair), mpl::net_screen_aznut(pair)};
        loom_check_bits(c, "scnfac.zhat", sf.zhat, lf.zhat);
        loom_check_bits(c, "scnfac.lzav", sf.lzav, lf.lzav);
        loom_check_bits(c, "scnfac.aznut", sf.aznut, lf.aznut);
        loom_check_bits(c, "screen_r", actual_screen(ps, sf), mpl::actual_screen_r(&pl, &lf, &rp));
        cmp_dual(c, "screen_d", actual_screen(psd, sf), mpl::actual_screen_d(&pld, &lf, &rp));
    }
}

// one-hot per species, all pairs 50/50, uniform, near-pure, and random
// (log-uniform over 10 decades) compositions, each summing to 1
inline std::vector<std::vector<Real>> compositions(std::mt19937_64& rng, int n_random)
{
    std::vector<std::vector<Real>> comps;
    for (int a = 0; a < NumSpec; ++a) {
        std::vector<Real> x(NumSpec, 0.0);
        x[a] = 1.0;
        comps.push_back(x);
        for (int b = a + 1; b < NumSpec; ++b) {
            std::vector<Real> y(NumSpec, 0.0);
            y[a] = 0.5;
            y[b] = 0.5;
            comps.push_back(y);
        }
    }
    comps.push_back(std::vector<Real>(NumSpec, 1.0 / NumSpec));
    {
        std::vector<Real> x(NumSpec, 1.0e-10);
        x[0] = 1.0 - (NumSpec - 1) * 1.0e-10;
        comps.push_back(x);
    }
    std::uniform_real_distribution<Real> u01(0.0, 1.0);
    for (int k = 0; k < n_random; ++k) {
        std::vector<Real> x(NumSpec);
        Real s = 0.0;
        for (auto& v : x) { v = std::pow(10.0, -10.0 * u01(rng)); s += v; }
        for (auto& v : x) { v /= s; }
        comps.push_back(x);
    }
    return comps;
}

} // namespace loomval_rates

inline long validate_rates()
{
    using namespace loomval_rates;
    LoomCheck c;
    c.module = "rates";

    std::vector<Real> temps;
    for (Real lt = 4.0; lt <= 16.0; lt += 0.01) { temps.push_back(std::pow(10.0, lt)); }
    std::vector<Real> branches;
    loomval_port::branch_temperatures(branches);
    for (Real tb0 : branches) {
        Real tb = tb0;
        for (int k = 0; k < 64; ++k) { tb = std::nextafter(tb, 0.0); }
        for (int k = 0; k < 129; ++k) { temps.push_back(tb); tb = std::nextafter(tb, 1.0e300); }
    }
    std::mt19937_64 rng(12345);
    std::uniform_real_distribution<Real> u(4.0, 16.0);
    for (int k = 0; k < 20000; ++k) { temps.push_back(std::pow(10.0, u(rng))); }

    loomval_port::save_params();
    for (int v = 0; v < loomval_port::num_param_variants; ++v) {
        loomval_port::apply_param_variant(v);
        for (Real T : temps) { one_T(c, T); }
    }
    loomval_port::restore_params();
    amrex::Print() << "LOOM_VALIDATE rates rate checks so far=" << c.checks << " (" << temps.size()
                   << " temperatures x " << loomval_port::num_param_variants << " parameter variants, "
                   << mpl::Rates::NumRates << " rates)\n";

    const auto comps = compositions(rng, 8);
    const bool save_skip = screening_rp::enable_debye_huckel_skip;
    for (int skip = 0; skip <= 1; ++skip) {
        screening_rp::enable_debye_huckel_skip = skip;
        for (const auto& X : comps) {
            for (Real lr = 2.0; lr <= 12.0; lr += 0.25) {
                for (Real lt = 4.0; lt <= 16.0; lt += 0.1) {
                    screen_one(c, std::pow(10.0, lt), std::pow(10.0, lr), X);
                }
            }
        }
    }
    screening_rp::enable_debye_huckel_skip = save_skip;

    const long m = loom_report(c);

    // sensitivity: T + 1 ulp on the Loom side must change some rate (not counted)
    {
        LoomCheck s;
        s.module = "rates-sensitivity";
        s.reported = 100;  // silence
        mpl::RatesParams rp{};
        fill_rates_params(&rp);
        for (Real T : {1.0e8, 3.0e8, 1.0e9, 5.0e9}) {
            DispatchAll<1>::run(s, T, std::nextafter(T, 1.0e300), &rp);
        }
        amrex::Print() << "LOOM_VALIDATE rates sensitivity (T +1ulp): " << s.mismatches << "/"
                       << s.checks << " differ (expected > 0)\n";
        if (s.mismatches == 0) {
            amrex::Print() << "LOOM_MISMATCH rates sensitivity: perturbation went undetected\n";
            return m + 1;
        }
    }
    return m;
}

#endif
