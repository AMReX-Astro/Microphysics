// Native bit-exact validation of dual.H against autodiff (host only).
// See README.md "Native validation" and the translation guide in dual.H.
//
// Part 1: every mpl dual primitive vs the autodiff expression it mirrors,
//         for autodiff::dual and autodiff::dual_array<1,3>, on edge values
//         (+-0, +-inf, NaN, denormals, +-DBL_MAX, ...) and random values over
//         wide magnitude ranges.
// Network-independent template: only the primitive tests. A port adds its
// own Part 2 (whole stock rate / screening / rhs.H expressions evaluated with
// autodiff vs the mpl call sequences the dual.H guide prescribes), which is
// network-specific and is not kept in the skill assets.
#ifndef VALIDATE_DUAL_H
#define VALIDATE_DUAL_H

#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <random>
#include <string>

#include <microphysics_autodiff.H>

#include <loom_validate.H>
#include "dual.H"

namespace loomval_dual {

using mpl::Real;
using mpl::Dual;
using mpl::Dual3;
using ad = autodiff::dual;
using ad3 = autodiff::dual_array<1, 3>;

// per-check mismatch tally (diagnostics)
inline std::map<std::string, long>& tally()
{
    static std::map<std::string, long> m;
    return m;
}

// Bit-exact compare, except that two NaNs whose sign/payload differ count as
// equal (and are tallied separately, see the report). LLVM leaves the sign
// and payload of a NaN result unspecified and freely rewrites e.g.
// (-x)*y -> -(x*y) in the inlined stock code, so NaN bits are a property of
// the optimizer, not of the source; any non-NaN difference is a mismatch.
inline long g_nan_payload_diffs = 0;

inline void cbits(LoomCheck& c, const char* what, double s, double m)
{
    if (std::isnan(s) && std::isnan(m)) {
        std::uint64_t a = 0;
        std::uint64_t b = 0;
        std::memcpy(&a, &s, sizeof a);
        std::memcpy(&b, &m, sizeof b);
        if (a != b) {
            ++g_nan_payload_diffs;
            tally()[std::string(what) + " [nan sign/payload only]"] += 1;
            m = s;
        }
    }
    if (!loom_check_bits(c, what, s, m)) {
        tally()[std::string(what) + " [value]"] += 1;
    }
}

// ---------------------------------------------------------------------------
// Ops: one interface over Dual / Dual3 so generic tests cover both
// ---------------------------------------------------------------------------

struct OpsD {
    using D = Dual;
    using AD = ad;
    static constexpr bool is3 = false;
    static D cst(Real n) { return mpl::dual_const(n); }
    static D seed(D d, int) { return d; }
    static D add(D a, D b) { return mpl::dual_add(a, b); }
    static D sub(D a, D b) { return mpl::dual_sub(a, b); }
    static D add_num(D a, Real n) { return mpl::dual_add_num(a, n); }
    static D add_ndm(D a, Real n, D b) { return mpl::dual_add_ndm(a, n, b); }
    static D mul(D a, D b) { return mpl::dual_mul(a, b); }
    static D mul_num(D a, Real n) { return mpl::dual_mul_num(a, n); }
    static D div(D a, D b) { return mpl::dual_div(a, b); }
    static D div_num(D a, Real n) { return mpl::dual_div_num(a, n); }
    static D neg(D a) { return mpl::dual_neg(a); }
    static D inv(D a) { return mpl::dual_inv(a); }
    static D exp(D a) { return mpl::dual_exp(a); }
    static D log(D a) { return mpl::dual_log(a); }
    static D log10(D a) { return mpl::dual_log10(a); }
    static D sqrt(D a) { return mpl::dual_sqrt(a); }
    static D cbrt(D a) { return mpl::dual_cbrt(a); }
    static D abs(D a) { return mpl::dual_abs(a); }
    static D pow_num(D a, Real n) { return mpl::dual_pow_num(a, n); }
    static D pow(D a, D b) { return mpl::dual_pow(a, b); }
    static D min(D a, D b) { return mpl::dual_min(a, b); }
    static D max(D a, D b) { return mpl::dual_max(a, b); }
    static D min_std(D a, D b) { return mpl::dual_min_std(a, b); }
    static D max_std(D a, D b) { return mpl::dual_max_std(a, b); }
    template <int N> static D powi(D a) { return mpl::dual_powi<N>(a); }

    static D make(Real v, const Real* g) { return D{v, g[0]}; }
    static AD to_ad(D d) { AD r; r.val = d.val; r.grad = d.grad; return r; }
    static void check(LoomCheck& c, const char* what, const AD& s, D m) {
        cbits(c, what, s.val, m.val);
        cbits(c, what, s.grad, m.grad);
    }
};

struct OpsD3 {
    using D = Dual3;
    using AD = ad3;
    static constexpr bool is3 = true;
    static D cst(Real n) { return mpl::dual3_const(n); }
    static D seed(D d, int k) { return mpl::dual3_seed(d, k); }
    static D add(D a, D b) { return mpl::dual3_add(a, b); }
    static D sub(D a, D b) { return mpl::dual3_sub(a, b); }
    static D add_num(D a, Real n) { return mpl::dual3_add_num(a, n); }
    static D add_ndm(D a, Real n, D b) { return mpl::dual3_add_ndm(a, n, b); }
    static D mul(D a, D b) { return mpl::dual3_mul(a, b); }
    static D mul_num(D a, Real n) { return mpl::dual3_mul_num(a, n); }
    static D div(D a, D b) { return mpl::dual3_div(a, b); }
    static D div_num(D a, Real n) { return mpl::dual3_div_num(a, n); }
    static D neg(D a) { return mpl::dual3_neg(a); }
    static D inv(D a) { return mpl::dual3_inv(a); }
    static D exp(D a) { return mpl::dual3_exp(a); }
    static D log(D a) { return mpl::dual3_log(a); }
    static D log10(D a) { return mpl::dual3_log10(a); }
    static D sqrt(D a) { return mpl::dual3_sqrt(a); }
    static D cbrt(D a) { return mpl::dual3_cbrt(a); }
    static D abs(D a) { return mpl::dual3_abs(a); }
    static D pow_num(D a, Real n) { return mpl::dual3_pow_num(a, n); }
    static D pow(D a, D b) { return mpl::dual3_pow(a, b); }
    static D min(D a, D b) { return mpl::dual3_min(a, b); }
    static D max(D a, D b) { return mpl::dual3_max(a, b); }
    static D min_std(D a, D b) { return mpl::dual3_min_std(a, b); }
    static D max_std(D a, D b) { return mpl::dual3_max_std(a, b); }
    template <int N> static D powi(D a) { return mpl::dual3_powi<N>(a); }

    static D make(Real v, const Real* g) {
        D r{};
        r.val = v;
        for (int i = 1; i <= 3; ++i) r.grad[i] = g[i - 1];
        return r;
    }
    static AD to_ad(const D& d) {
        AD r;
        r.val = d.val;
        for (int i = 1; i <= 3; ++i) r.grad(i) = d.grad[i];
        return r;
    }
    static void check(LoomCheck& c, const char* what, const AD& s, const D& m) {
        cbits(c, what, s.val, m.val);
        for (int i = 1; i <= 3; ++i) cbits(c, what, s.grad(i), m.grad[i]);
        cbits(c, what, 0.0, m.grad[0]);
    }
};

// ---------------------------------------------------------------------------
// inputs
// ---------------------------------------------------------------------------

inline const std::vector<Real>& edge_values()
{
    static const std::vector<Real> v = {
        0.0, -0.0, 1.0, -1.0, 0.5, -2.0, 3.0, 1.0e-310, -1.0e-310,
        1.0e-300, -1.0e300, std::numeric_limits<Real>::max(),
        -std::numeric_limits<Real>::max(), std::numeric_limits<Real>::infinity(),
        -std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN(),
        0.1, 7.25e-5, -1.6
    };
    return v;
}

struct Rng {
    std::mt19937_64 gen{0x5eed1234abcdULL};
    Real uni(Real a, Real b) { return std::uniform_real_distribution<Real>(a, b)(gen); }
    int pick(int n) { return std::uniform_int_distribution<int>(0, n - 1)(gen); }
    // sign * mantissa * 10^e over a wide or a moderate exponent range
    Real wide() {
        const int mode = pick(10);
        if (mode == 0) return edge_values()[pick(static_cast<int>(edge_values().size()))];
        const Real e = (mode < 4) ? uni(-300.0, 300.0) : uni(-12.0, 12.0);
        const Real s = (pick(3) == 0) ? -1.0 : 1.0;
        return s * uni(1.0, 10.0) * std::pow(10.0, e);
    }
    Real logu(Real lo, Real hi) { return std::pow(10.0, uni(std::log10(lo), std::log10(hi))); }
};

// ---------------------------------------------------------------------------
// Part 1: primitives
// ---------------------------------------------------------------------------

template <class O>
void prim_tests(LoomCheck& c, const typename O::D& a, const typename O::D& b, Real n)
{
    using AD = typename O::AD;
    const AD A = O::to_ad(a);
    const AD B = O::to_ad(b);

    { AD r = n; O::check(c, "const", r, O::cst(n)); }
    { AD r = A + B; O::check(c, "add", r, O::add(a, b)); }
    { AD r = A - B; O::check(c, "sub", r, O::sub(a, b)); }
    { AD r = A + n; O::check(c, "add_num", r, O::add_num(a, n)); }
    { AD r = A - n; O::check(c, "sub_num->add_num(-n)", r, O::add_num(a, -n)); }
    { AD r = n + A; O::check(c, "num+dual", r, O::add(O::cst(n), a)); }
    { AD r = n - A; O::check(c, "num-dual", r, O::sub(O::cst(n), a)); }
    { AD r = A + n * B; O::check(c, "add_ndm", r, O::add_ndm(a, n, b)); }
    { AD r = A - n * B; O::check(c, "sub_ndm->add_ndm(-n)", r, O::add_ndm(a, -n, b)); }
    { AD r = A + B * n; O::check(c, "dual+dual*n", r, O::add(a, O::mul_num(b, n))); }
    { AD r = A * B; O::check(c, "mul", r, O::mul(a, b)); }
    { AD r = A * n; O::check(c, "mul_num", r, O::mul_num(a, n)); }
    { AD r = n * A; O::check(c, "ndm", r, O::mul_num(a, n)); }
    { AD r = n * A * B; O::check(c, "n*a*b", r, O::mul(O::mul_num(a, n), b)); }
    { AD r = A * (B * B); O::check(c, "a*(b*b)", r, O::mul(a, O::mul(b, b))); }
    { AD r = A * (n * B); O::check(c, "a*(n*b)", r, O::mul(a, O::mul_num(b, n))); }
    { AD r = A / B; O::check(c, "div", r, O::div(a, b)); }
    { AD r = A / n; O::check(c, "div_num", r, O::div_num(a, n)); }
    { AD r = n / A; O::check(c, "num/dual", r, O::div(O::cst(n), a)); }
    { AD r = A / (n * B); O::check(c, "a/(n*b)", r, O::div(a, O::mul_num(b, n))); }
    { AD r = A / (A + B); O::check(c, "a/(a+b)", r, O::div(a, O::add(a, b))); }
    { AD r = -A; O::check(c, "neg", r, O::neg(a)); }
    { AD r = -(n * A); O::check(c, "-(n*a)", r, O::mul_num(a, -n)); }
    { AD r = n * (-A); O::check(c, "n*(-a)", r, O::mul_num(a, -n)); }
    { AD r = (-A) * (-B); O::check(c, "(-a)*(-b)", r, O::mul(a, b)); }
    { AD r = A * (-B); O::check(c, "a*(-b)", r, O::neg(O::mul(a, b))); }
    { AD r = A / (-B); O::check(c, "a/(-b)", r, O::neg(O::div(a, b))); }
    { AD r = (-A) + (-B); O::check(c, "(-a)+(-b)", r, O::neg(O::add(a, b))); }
    { AD r = -A + B; O::check(c, "-a+b", r, O::add(O::neg(a), b)); }
    { AD r = A - (-B); O::check(c, "a-(-b)", r, O::add(a, b)); }
    { AD r = autodiff::detail::inverse(A); O::check(c, "inv", r, O::inv(a)); }
    { AD r = admath::exp(A); O::check(c, "exp", r, O::exp(a)); }
    { AD r = admath::log(A); O::check(c, "log", r, O::log(a)); }
    { AD r = admath::log10(A); O::check(c, "log10", r, O::log10(a)); }
    { AD r = admath::sqrt(A); O::check(c, "sqrt", r, O::sqrt(a)); }
    { AD r = admath::cbrt(A); O::check(c, "cbrt", r, O::cbrt(a)); }
    { AD r = admath::abs(A); O::check(c, "abs", r, O::abs(a)); }
    { AD r = admath::pow(A, n); O::check(c, "pow_num", r, O::pow_num(a, n)); }
    { AD r = admath::pow(A, B); O::check(c, "pow", r, O::pow(a, b)); }
    { AD r = admath::pow(n, B); O::check(c, "num^dual", r, O::pow(O::cst(n), b)); }
    { AD r = n * admath::exp(A); O::check(c, "n*exp(a)", r, O::mul(O::cst(n), O::exp(a))); }
    { AD r = admath::exp(-admath::pow(n * A, 3.263)); O::check(c, "exp(-pow(n*a))", r, O::exp(O::neg(O::pow_num(O::mul_num(a, n), 3.263)))); }
    { AD r = admath::sqrt(admath::sqrt(A)); O::check(c, "sqrt(sqrt)", r, O::sqrt(O::sqrt(a))); }
    { AD r = admath::powi<-3>(A); O::check(c, "powi-3", r, O::template powi<-3>(a)); }
    { AD r = admath::powi<-1>(A); O::check(c, "powi-1", r, O::template powi<-1>(a)); }
    { AD r = admath::powi<0>(A); O::check(c, "powi0", r, O::template powi<0>(a)); }
    { AD r = admath::powi<1>(A); O::check(c, "powi1", r, O::template powi<1>(a)); }
    { AD r = admath::powi<2>(A); O::check(c, "powi2", r, O::template powi<2>(a)); }
    { AD r = admath::powi<3>(A); O::check(c, "powi3", r, O::template powi<3>(a)); }
    { AD r = admath::powi<4>(A); O::check(c, "powi4", r, O::template powi<4>(a)); }
    { AD r = admath::powi<6>(A); O::check(c, "powi6", r, O::template powi<6>(a)); }
    { AD r = admath::powi<7>(A); O::check(c, "powi7", r, O::template powi<7>(a)); }
    // const lvalues: amrex::min/max (std rule); non-const or expressions: autodiff
    { AD r = admath::min(A, B); O::check(c, "min(const,const)", r, O::min_std(a, b)); }
    { AD r = admath::max(A, B); O::check(c, "max(const,const)", r, O::max_std(a, b)); }
    { AD x = A; AD y = B; AD r = admath::min(x, y); O::check(c, "min", r, O::min(a, b)); }
    { AD x = A; AD y = B; AD r = admath::max(x, y); O::check(c, "max", r, O::max(a, b)); }
    { AD r = admath::min(A + 0.0, B); O::check(c, "min(expr,dual)", r, O::min(O::add_num(a, 0.0), b)); }
    { AD r = admath::max(A * B, B); O::check(c, "max(expr,dual)", r, O::max(O::mul(a, b), b)); }
    { AD r = admath::min(n, A); O::check(c, "min(n,a)", r, O::min(O::cst(n), a)); }
    { AD r = admath::max(A, n); O::check(c, "max(a,n)", r, O::max(a, O::cst(n))); }
    { AD r = A; r += B; O::check(c, "+=dual", r, O::add(a, b)); }
    { AD r = A; r -= B; O::check(c, "-=dual", r, O::sub(a, b)); }
    { AD r = A; r *= B; O::check(c, "*=dual", r, O::mul(a, b)); }
    { AD r = A; r /= B; O::check(c, "/=dual", r, O::div(a, b)); }
    { AD r = A; r += n; O::check(c, "+=num", r, O::add(a, O::cst(n))); }
    { AD r = A; r -= n; O::check(c, "-=num", r, O::sub(a, O::cst(n))); }
    { AD r = A; r *= n; O::check(c, "*=num", r, O::mul(a, O::cst(n))); }
    { AD r = A; r /= n; O::check(c, "/=num", r, O::div(a, O::cst(n))); }
    { AD r = A; r *= A; O::check(c, "*=self", r, O::mul(a, a)); }
    { AD r = A; r += B * n; O::check(c, "+=expr", r, O::add(a, O::mul_num(b, n))); }
    { AD r = A; r = r * B + n; O::check(c, "r=r*b+n", r, O::add_num(O::mul(a, b), n)); }
    loom_check_int(c, "cmp>", (A > n) ? 1 : 0, (a.val > n) ? 1 : 0);
    loom_check_int(c, "cmp<=", (A <= B) ? 1 : 0, (a.val <= b.val) ? 1 : 0);
    loom_check_int(c, "cmp>=", (A >= n) ? 1 : 0, (a.val >= n) ? 1 : 0);
}

template <class O>
void run_prims(LoomCheck& c, Rng& rng)
{
    using D = typename O::D;
    const auto& ev = edge_values();
    const int ne = static_cast<int>(ev.size());
    Real g[3] = {0.0, 0.0, 0.0};
    // all (val_a, val_b, n) edge triples, with edge gradients
    for (int i = 0; i < ne; ++i) {
        for (int j = 0; j < ne; ++j) {
            for (int k = 0; k < ne; ++k) {
                for (int m = 0; m < 3; ++m) g[m] = ev[(i + j + k + 5 * m) % ne];
                const D a = O::make(ev[i], g);
                for (int m = 0; m < 3; ++m) g[m] = ev[(3 * i + k + 7 * m + 1) % ne];
                const D b = O::make(ev[j], g);
                prim_tests<O>(c, a, b, ev[k]);
            }
        }
    }
    // random over wide and moderate magnitudes
    for (int t = 0; t < 40000; ++t) {
        for (int m = 0; m < 3; ++m) g[m] = (rng.pick(8) == 0) ? 0.0 : rng.wide();
        const D a = O::make(rng.wide(), g);
        for (int m = 0; m < 3; ++m) g[m] = (rng.pick(8) == 0) ? -0.0 : rng.wide();
        const D b = O::make(rng.wide(), g);
        prim_tests<O>(c, a, b, rng.wide());
    }
    // seeds
    if constexpr (O::is3) {
        for (int k = 1; k <= 3; ++k) {
            ad3 x = 2.5;
            ad3 y = 3.5;
            ad3 z = 4.5;
            autodiff::seed_array(x, y, z);
            const Dual3 s = mpl::dual3_seed(mpl::dual3_const(1.5 + k), k);
            ad3 w = (k == 1) ? x : ((k == 2) ? y : z);
            O::check(c, "seed_array", w, s);
        }
    } else {
        ad x = 2.5;
        autodiff::seed(x);
        O::check(c, "seed", x, mpl::dual_seed(2.5));
    }
}

} // namespace loomval_dual

inline long validate_dual()
{
    using namespace loomval_dual;
    LoomCheck c;
    c.module = "dual";
    Rng rng;
    run_prims<OpsD>(c, rng);
    run_prims<OpsD3>(c, rng);
    amrex::Print() << "LOOM_VALIDATE dual NaN-vs-NaN sign/payload-only differences (not counted): "
                   << g_nan_payload_diffs << "\n";
    for (const auto& kv : tally()) {
        amrex::Print() << "LOOM_VALIDATE dual mismatch tally: " << kv.first << " = " << kv.second << "\n";
    }
    return loom_report(c);
}

#endif
