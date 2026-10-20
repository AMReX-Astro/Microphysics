// Native bit-exact validation of the translated helmholtz EOS (eos.H)
// against the stock Microphysics eos() and its helpers. Host only; included
// by loom_validate.cpp under -DLOOM_VALIDATE_EOS. See README.md.
#ifndef VALIDATE_EOS_H
#define VALIDATE_EOS_H

#include <cmath>
#include <cstdint>
#include <cstring>
#include <random>
#include <string>
#include <type_traits>
#include <vector>

#include <extern_parameters.H>
#include <eos.H>            // stock (interfaces/eos.H; the kernel dir is last on -I)
#include <react_util.H>
#include <loom_validate.H>

#include "eos.H"            // translated (this directory)
#include <loom_eos_params.H>   // host fill_eos_* (unit_test/test_react); not loom_params.H,
                               // so this builds before rates.H / network.H exist

static_assert(std::is_same_v<amrex::Real, double>, "Loom kernel needs PRECISION=DOUBLE");
static_assert(NumSpec == mpl::eos_nspec, "mpl::eos_nspec must match NumSpec");
#if NAUX_NET > 0
#error "mpl EOS translation assumes NAUX_NET = 0"
#endif
#ifdef AUX_THERMO
#error "mpl EOS translation assumes no AUX_THERMO"
#endif

// The translation resolves every has_*<T> for this trait set; the burn's
// two EOS state types must have it.
template <typename T>
constexpr bool loom_eos_trait_set_ok =
    has_xn<T> && has_energy<T> && has_eta<T> && has_dedA<T> && has_dedZ<T> &&
    !has_pressure<T> && !has_entropy<T> && !has_enthalpy<T> &&
    !has_dpdA<T> && !has_dpdZ<T> && !has_xne_xnp<T> && !has_pele_ppos<T>;
static_assert(loom_eos_trait_set_ok<burn_t>, "burn_t EOS trait set changed");
static_assert(loom_eos_trait_set_ok<eos_re_extra_t>, "eos_re_extra_t EOS trait set changed");

// stock helmholtz constants used by the translation
static_assert(helmholtz::imax == mpl::helmholtz::imax && helmholtz::jmax == mpl::helmholtz::jmax);
static_assert(helmholtz::tstpi == mpl::helmholtz::tstpi && helmholtz::dstpi == mpl::helmholtz::dstpi);
static_assert(helmholtz::tlo == mpl::helmholtz::tlo && helmholtz::dlo == mpl::helmholtz::dlo);
static_assert(helmholtz::avo_eos == mpl::helmholtz::avo_eos && helmholtz::kerg == mpl::helmholtz::kerg);
static_assert(std::numbers::pi_v<double> == mpl::eos_pi);
static_assert(static_cast<int>(eos_input_rt) == mpl::eos_input_rt &&
              static_cast<int>(eos_input_rh) == mpl::eos_input_rh &&
              static_cast<int>(eos_input_tp) == mpl::eos_input_tp &&
              static_cast<int>(eos_input_rp) == mpl::eos_input_rp &&
              static_cast<int>(eos_input_re) == mpl::eos_input_re &&
              static_cast<int>(eos_input_ps) == mpl::eos_input_ps &&
              static_cast<int>(eos_input_ph) == mpl::eos_input_ph &&
              static_cast<int>(eos_input_th) == mpl::eos_input_th);
static_assert(static_cast<int>(EOS::itemp) == mpl::EOS::itemp &&
              static_cast<int>(EOS::idens) == mpl::EOS::idens &&
              static_cast<int>(EOS::iener) == mpl::EOS::iener &&
              static_cast<int>(EOS::ienth) == mpl::EOS::ienth &&
              static_cast<int>(EOS::ientr) == mpl::EOS::ientr &&
              static_cast<int>(EOS::ipres) == mpl::EOS::ipres);

// fill_eos_params / fill_eos_tables: <loom_eos_params.H>

namespace loom_eos_val {

using mpl::Real;

// The EOS-visible fields of burn_t (the kernel's burn_t mirror must have
// them). Lives in its own namespace so ADL never pulls in the stock
// global-namespace eos templates.
struct BurnEos {
    Real rho{};
    Real T{};
    Real e{};
    Real xn[NumSpec]{};
    Real dedr{};
    Real dedT{};
    Real dedA{};
    Real dedZ{};
    Real cv{};
    Real mu{};
    Real mu_e{};
    Real y_e{};
    Real eta{};
    Real abar{};
    Real zbar{};
};

struct Ctx {
    std::string tag;
    int input{-1};
    Real rho{}, T{}, e{};
};

inline bool same_bits(Real a, Real b)
{
    std::uint64_t x = 0, y = 0;
    std::memcpy(&x, &a, sizeof x);
    std::memcpy(&y, &b, sizeof y);
    return x == y;
}

inline void chk(LoomCheck& c, const Ctx& ctx, const char* field, Real stock, Real loom)
{
    if (same_bits(stock, loom)) {
        ++c.checks;
        return;
    }
    char buf[256];
    std::snprintf(buf, sizeof buf, "%s:%s input=%d rho=%.17g T=%.17g e=%.17g",
                  ctx.tag.c_str(), field, ctx.input, ctx.rho, ctx.T, ctx.e);
    loom_check_bits(c, buf, stock, loom);
}

template <typename S, typename M>
void copy_in(const S& s, M* m)
{
    m->rho = s.rho;
    m->T = s.T;
    m->e = s.e;
    for (int n = 0; n < NumSpec; ++n) {
        m->xn[n] = s.xn[n];
    }
    m->dedr = s.dedr;
    m->dedT = s.dedT;
    m->dedA = s.dedA;
    m->dedZ = s.dedZ;
    m->cv = s.cv;
    m->mu = s.mu;
    m->mu_e = s.mu_e;
    m->y_e = s.y_e;
    m->eta = s.eta;
    m->abar = s.abar;
    m->zbar = s.zbar;
}

// every field of the trait set (all of eos_re_extra_t, all EOS-visible
// fields of burn_t)
template <typename S, typename M>
void check_all(LoomCheck& c, const Ctx& ctx, const S& s, const M& m)
{
    chk(c, ctx, "rho", s.rho, m.rho);
    chk(c, ctx, "T", s.T, m.T);
    chk(c, ctx, "e", s.e, m.e);
    for (int n = 0; n < NumSpec; ++n) {
        chk(c, ctx, "xn", s.xn[n], m.xn[n]);
    }
    chk(c, ctx, "dedr", s.dedr, m.dedr);
    chk(c, ctx, "dedT", s.dedT, m.dedT);
    chk(c, ctx, "dedA", s.dedA, m.dedA);
    chk(c, ctx, "dedZ", s.dedZ, m.dedZ);
    chk(c, ctx, "cv", s.cv, m.cv);
    chk(c, ctx, "mu", s.mu, m.mu);
    chk(c, ctx, "mu_e", s.mu_e, m.mu_e);
    chk(c, ctx, "y_e", s.y_e, m.y_e);
    chk(c, ctx, "eta", s.eta, m.eta);
    chk(c, ctx, "abar", s.abar, m.abar);
    chk(c, ctx, "zbar", s.zbar, m.zbar);
}

// Fill a stock state with (rho, T, e, xn) and arbitrary values in every
// other field (the EOS must overwrite or ignore them identically).
template <typename S>
S make_state(Real rho, Real T, Real e, const Real* xn, Real junk)
{
    S s{};
    s.rho = rho;
    s.T = T;
    s.e = e;
    for (int n = 0; n < NumSpec; ++n) {
        s.xn[n] = xn[n];
    }
    s.dedr = junk;
    s.dedT = 2.0 * junk;
    s.dedA = 3.0 * junk;
    s.dedZ = 4.0 * junk;
    s.cv = 5.0 * junk;
    s.mu = 6.0 * junk;
    s.mu_e = 7.0 * junk;
    s.y_e = 8.0 * junk;
    s.eta = 9.0 * junk;
    s.abar = 10.0 * junk;
    s.zbar = 11.0 * junk;
    return s;
}

// stock eos() vs mpl::eos() on the same input, all fields compared
template <typename S, typename M>
void run_eos(LoomCheck& c, const char* tag, int input, const S& init,
             const mpl::EosParams* p, const mpl::EosTables* t)
{
    S s = init;
    M m{};
    copy_in(init, &m);
    eos(static_cast<eos_input_t>(input), s);
    mpl::eos(input, &m, *p, *t);
    Ctx ctx{tag, input, init.rho, init.T, init.e};
    check_all(c, ctx, s, m);
}

// stock reset_inputs vs mpl::reset_inputs (after composition), including
// the has_been_reset flag
template <typename S, typename M>
void run_reset(LoomCheck& c, const char* tag, int input, const S& init,
               const mpl::EosParams* p, const mpl::EosTables* t)
{
    S s = init;
    M m{};
    copy_in(init, &m);
    composition(s);
    mpl::composition(&m);
    bool hs = false;
    bool hm = false;
    reset_inputs(static_cast<eos_input_t>(input), s, hs);
    mpl::reset_inputs(input, &m, hm, *p, *t);
    Ctx ctx{tag, input, init.rho, init.T, init.e};
    check_all(c, ctx, s, m);
    loom_check_int(c, (std::string(tag) + ":has_been_reset").c_str(), hs ? 1 : 0, hm ? 1 : 0);
}

// the physics pieces one at a time, in actual_eos order
inline void run_pieces(LoomCheck& c, const Real* xn, Real rho, Real T, Real junk,
                       const mpl::EosParams* p, const mpl::EosTables* t)
{
    burn_t s = make_state<burn_t>(rho, T, 0.0, xn, junk);
    BurnEos m{};
    copy_in(s, &m);
    Ctx ctx{"pieces", -1, rho, T, 0.0};

    composition(s);
    mpl::composition(&m);
    ctx.tag = "composition";
    check_all(c, ctx, s, m);

    apply_radiation(s);
    mpl::apply_radiation(&m, *p);
    ctx.tag = "apply_radiation";
    check_all(c, ctx, s, m);

    apply_ions(s);
    mpl::apply_ions(&m);
    ctx.tag = "apply_ions";
    check_all(c, ctx, s, m);

    apply_electrons(s);
    mpl::apply_electrons(&m, *t);
    ctx.tag = "apply_electrons";
    check_all(c, ctx, s, m);

    apply_coulomb_corrections(s);
    mpl::apply_coulomb_corrections(&m);
    ctx.tag = "apply_coulomb_corrections";
    check_all(c, ctx, s, m);

    // coulomb on a state with an arbitrary (possibly negative) energy,
    // to exercise the "disable if e goes negative" branch
    s.e = -junk;
    m.e = -junk;
    apply_coulomb_corrections(s);
    mpl::apply_coulomb_corrections(&m);
    ctx.tag = "apply_coulomb_corrections(e<0)";
    check_all(c, ctx, s, m);

    // finalize_state for the re mode
    finalize_state(eos_input_re, s, 1.25 * junk, 0.0, 0.0);
    mpl::finalize_state(mpl::eos_input_re, &m, 1.25 * junk, 0.0, 0.0, *p);
    ctx.tag = "finalize_state";
    check_all(c, ctx, s, m);
}

// composition_derivatives on eos_re_extra_t after eos(eos_input_re), as in
// the Strang jac (correct_jacobian_for_const_e)
inline void run_xderivs(LoomCheck& c, const eos_re_extra_t& init,
                        const mpl::EosParams* p, const mpl::EosTables* t)
{
    eos_re_extra_t s = init;
    mpl::eos_re_extra_t m{};
    copy_in(init, &m);
    eos(eos_input_re, s);
    mpl::eos(mpl::eos_input_re, &m, *p, *t);
    eos_xderivs_t xs = composition_derivatives(s);
    mpl::eos_xderivs_t xm = mpl::composition_derivatives(m);
    Ctx ctx{"composition_derivatives", mpl::eos_input_re, init.rho, init.T, init.e};
    for (int n = 0; n < NumSpec; ++n) {
        chk(c, ctx, "dedX", xs.dedX[n], xm.dedX[n]);
        chk(c, ctx, "dpdX", xs.dpdX[n], xm.dpdX[n]);
        chk(c, ctx, "dhdX", xs.dhdX[n], xm.dhdX[n]);
    }
}

inline void run_polys(LoomCheck& c, std::mt19937_64& rng)
{
    Ctx ctx{"psi", -1, 0.0, 0.0, 0.0};
    std::uniform_real_distribution<Real> u(-0.5, 1.5);
    for (int k = 0; k < 20001; ++k) {
        Real z = (k < 10001) ? (-0.5 + 2.0 * k / 10000.0) : u(rng);
        ctx.rho = z;
        chk(c, ctx, "psi0", psi0(z), mpl::psi0(z));
        chk(c, ctx, "dpsi0", dpsi0(z), mpl::dpsi0(z));
        chk(c, ctx, "ddpsi0", ddpsi0(z), mpl::ddpsi0(z));
        chk(c, ctx, "psi1", psi1(z), mpl::psi1(z));
        chk(c, ctx, "dpsi1", dpsi1(z), mpl::dpsi1(z));
        chk(c, ctx, "ddpsi1", ddpsi1(z), mpl::ddpsi1(z));
        chk(c, ctx, "psi2", psi2(z), mpl::psi2(z));
        chk(c, ctx, "dpsi2", dpsi2(z), mpl::dpsi2(z));
        chk(c, ctx, "ddpsi2", ddpsi2(z), mpl::ddpsi2(z));
        chk(c, ctx, "xpsi0", xpsi0(z), mpl::xpsi0(z));
        chk(c, ctx, "xpsi1", xpsi1(z), mpl::xpsi1(z));
    }
    std::uniform_real_distribution<Real> v(-1.0e3, 1.0e3);
    for (int k = 0; k < 2000; ++k) {
        Real fi[36];
        Real wt[6];
        for (Real& x : fi) x = v(rng);
        for (Real& x : wt) x = v(rng);
        Real a[6];
        Real b[6];
        fwt(fi, wt, a);
        mpl::fwt(fi, wt, b);
        for (int i = 0; i < 6; ++i) {
            chk(c, ctx, "fwt", a[i], b[i]);
        }
    }
}

// One table lookup by hand: the flat-index macros vs the stock arrays.
inline void run_table_index(LoomCheck& c, const mpl::EosTables* t)
{
    Ctx ctx{"table_index", -1, 0.0, 0.0, 0.0};
    for (int j = 0; j < helmholtz::jmax; j += 7) {
        for (int i = 0; i < helmholtz::imax; i += 5) {
            for (int m = 0; m < 9; ++m) {
                chk(c, ctx, "f", helmholtz::f[j][i][m], t->f[MPL_HELM_F(j, i, m)]);
            }
            for (int m = 0; m < 4; ++m) {
                chk(c, ctx, "ef", helmholtz::ef[j][i][m], t->ef[MPL_HELM_F4(j, i, m)]);
                chk(c, ctx, "dpdf", helmholtz::dpdf[j][i][m], t->dpdf[MPL_HELM_F4(j, i, m)]);
                chk(c, ctx, "xf", helmholtz::xf[j][i][m], t->xf[MPL_HELM_F4(j, i, m)]);
            }
        }
    }
}

// The compositions test_react uses (setup_composition / get_xn for
// n_cell zones, and the uniform one), single-species ones, and random ones.
inline std::vector<std::vector<Real>> make_compositions(std::mt19937_64& rng, int n_random,
                                                        int* n_test_react)
{
    std::vector<std::vector<Real>> comps;
    const int nz = 16;   // test_react n_cell (16 in the stock inputs files)
    init_t cd = setup_composition(nz);
    for (int k = 0; k < nz; ++k) {
        std::vector<Real> xn(NumSpec);
        get_xn(k, cd, xn.data(), false);
        comps.push_back(xn);
    }
    {
        std::vector<Real> xn(NumSpec);
        get_xn(0, cd, xn.data(), true);
        comps.push_back(xn);
    }
    *n_test_react = static_cast<int>(comps.size());

    for (int s = 0; s < NumSpec; ++s) {
        std::vector<Real> xn(NumSpec, 0.0);
        xn[s] = 1.0;
        comps.push_back(xn);
    }
    std::uniform_real_distribution<Real> u(0.0, 1.0);
    for (int r = 0; r < n_random; ++r) {
        std::vector<Real> xn(NumSpec);
        Real sum = 0.0;
        for (int n = 0; n < NumSpec; ++n) {
            Real w = u(rng);
            // sparse and very small abundances too
            if (r % 3 == 1 && u(rng) < 0.5) w = 0.0;
            if (r % 3 == 2) w = std::pow(w, 8.0);
            xn[n] = w;
            sum += w;
        }
        if (sum == 0.0) {
            xn[0] = 1.0;
            sum = 1.0;
        }
        for (int n = 0; n < NumSpec; ++n) xn[n] /= sum;
        comps.push_back(xn);
    }
    return comps;
}

inline std::vector<Real> log_grid(Real lo, Real hi, Real step)
{
    std::vector<Real> v;
    int n = static_cast<int>(std::floor((hi - lo) / step + 0.5));
    for (int i = 0; i <= n; ++i) {
        v.push_back(std::pow(10.0, lo + i * step));
    }
    return v;
}

// rho and T samples: a coarse grid well beyond the table / reset limits,
// a dense grid over the test_react range, and the exact limits.
inline void make_grids(std::vector<Real>& rhos, std::vector<Real>& temps, bool dense, bool sparse = false)
{
    rhos = log_grid(-14.0, 16.5, sparse ? 1.0 : 0.5);
    temps = log_grid(1.0, 16.5, sparse ? 0.5 : 0.25);
    if (dense) {
        auto r2 = log_grid(5.5, 9.5, 0.05);     // rho 1e6 - 1e9 and margin
        auto t2 = log_grid(5.5, 10.5, 0.05);
        rhos.insert(rhos.end(), r2.begin(), r2.end());
        temps.insert(temps.end(), t2.begin(), t2.end());
        for (Real x : {1.0e9 * 3.0, 1.0e10, 1.0e11, 1.0e12}) rhos.push_back(x);
        for (Real x : {1.0e11, 1.0e12, 1.0e13, 1.0e14, 1.0e15}) temps.push_back(x);
        // T on exact table nodes (xt = 0) and just below the upper limit
        for (int j = 0; j < helmholtz::jmax; j += 20) temps.push_back(helmholtz::t[j]);
    }
    for (Real x : {EOSData::mindens, EOSData::maxdens,
                   std::nextafter(EOSData::mindens, 0.0), std::nextafter(EOSData::maxdens, 1.0e300),
                   0.5 * EOSData::mindens, 2.0 * EOSData::maxdens}) {
        rhos.push_back(x);
    }
    for (Real x : {EOSData::mintemp, EOSData::maxtemp,
                   std::nextafter(EOSData::mintemp, 0.0), std::nextafter(EOSData::maxtemp, 1.0e300),
                   0.5 * EOSData::mintemp, 2.0 * EOSData::maxtemp}) {
        temps.push_back(x);
    }
}

struct Counts {
    long rt{0}, re{0}, other{0}, reset{0}, pieces{0}, xderivs{0};
};

// The main sweep for one parameter set.
inline void sweep(LoomCheck& c, const char* tag, bool dense, int comp_stride,
                  const std::vector<std::vector<Real>>& comps, int n_test_react,
                  std::mt19937_64& rng, Counts& cnt)
{
    mpl::EosParams p{};
    fill_eos_params(&p);
    mpl::EosTables t{};
    fill_eos_tables(&t);

    std::vector<Real> rhos;
    std::vector<Real> temps;
    std::uniform_real_distribution<Real> uf(-0.6, 0.6);
    std::uniform_real_distribution<Real> uj(-1.0e3, 1.0e3);

    for (int ic = 0; ic < static_cast<int>(comps.size()); ic += comp_stride) {
        const Real* xn = comps[ic].data();
        // test_react compositions on the dense grid, the others coarse
        make_grids(rhos, temps, dense && ic < n_test_react);
        for (Real rho : rhos) {
            for (Real T : temps) {
                Real junk = uj(rng);

                // eos_input_rt (integrator_setup), both state types
                run_eos<burn_t, BurnEos>(c, tag, mpl::eos_input_rt,
                                         make_state<burn_t>(rho, T, junk, xn, junk), &p, &t);
                run_eos<eos_re_extra_t, mpl::eos_re_extra_t>(c, tag, mpl::eos_input_rt,
                                         make_state<eos_re_extra_t>(rho, T, junk, xn, junk), &p, &t);
                cnt.rt += 2;

                // eos_input_re (update_thermodynamics, jac): target energy
                // from a stock rt call at a different temperature, so the
                // Newton iteration runs from T
                Real T2 = T * std::pow(10.0, uf(rng));
                burn_t tgt = make_state<burn_t>(rho, T2, 0.0, xn, 0.0);
                eos(eos_input_rt, tgt);
                Real e_target = tgt.e;
                run_eos<burn_t, BurnEos>(c, tag, mpl::eos_input_re,
                                         make_state<burn_t>(rho, T, e_target, xn, junk), &p, &t);
                run_eos<eos_re_extra_t, mpl::eos_re_extra_t>(c, tag, mpl::eos_input_re,
                                         make_state<eos_re_extra_t>(rho, T, e_target, xn, junk), &p, &t);
                cnt.re += 2;

                // individual pieces and composition_derivatives
                run_pieces(c, xn, rho, T, junk, &p, &t);
                cnt.pieces += 1;
                run_xderivs(c, make_state<eos_re_extra_t>(rho, T, e_target, xn, junk), &p, &t);
                cnt.xderivs += 1;
            }
        }
    }
}

// Edge cases reset_inputs / eos_reset / the iteration branch on, every
// input mode (including the ones whose bodies are empty for this trait
// set), both state types.
inline void edges(LoomCheck& c, const char* tag, const std::vector<std::vector<Real>>& comps,
                  std::mt19937_64& rng, Counts& cnt)
{
    mpl::EosParams p{};
    fill_eos_params(&p);
    mpl::EosTables t{};
    fill_eos_tables(&t);

    std::vector<Real> rhos;
    std::vector<Real> temps;
    make_grids(rhos, temps, false, true);
    std::uniform_real_distribution<Real> uj(-1.0e3, 1.0e3);

    const Real e_edges[] = {-1.0e18, -1.0, -0.0, 0.0, 0.5 * EOSData::mine, EOSData::mine,
                            std::nextafter(EOSData::mine, 0.0), 1.0e10, 1.0e17, 1.0e25,
                            EOSData::maxe, std::nextafter(EOSData::maxe, 1.0e300), 2.0 * EOSData::maxe,
                            std::numeric_limits<Real>::infinity()};

    for (int ic = 0; ic < static_cast<int>(comps.size()); ic += 6) {
        const Real* xn = comps[ic].data();
        for (Real rho : rhos) {
            for (Real T : temps) {
                Real junk = uj(rng);
                for (Real e : e_edges) {
                    for (int input = 0; input <= 7; ++input) {
                        run_reset<burn_t, BurnEos>(c, tag, input,
                                                   make_state<burn_t>(rho, T, e, xn, junk), &p, &t);
                        run_reset<eos_re_extra_t, mpl::eos_re_extra_t>(c, tag, input,
                                                   make_state<eos_re_extra_t>(rho, T, e, xn, junk), &p, &t);
                        cnt.reset += 2;
                    }
                    // the burn's modes on edge energies
                    run_eos<burn_t, BurnEos>(c, tag, mpl::eos_input_re,
                                             make_state<burn_t>(rho, T, e, xn, junk), &p, &t);
                    run_eos<eos_re_extra_t, mpl::eos_re_extra_t>(c, tag, mpl::eos_input_re,
                                             make_state<eos_re_extra_t>(rho, T, e, xn, junk), &p, &t);
                    cnt.re += 2;
                }
                // modes the burn does not use (empty branches for this trait
                // set, but statically reachable through the runtime input)
                for (int input : {mpl::eos_input_rh, mpl::eos_input_tp, mpl::eos_input_rp,
                                  mpl::eos_input_ps, mpl::eos_input_ph, mpl::eos_input_th}) {
                    run_eos<burn_t, BurnEos>(c, tag, input,
                                             make_state<burn_t>(rho, T, 1.0e17, xn, junk), &p, &t);
                    run_eos<eos_re_extra_t, mpl::eos_re_extra_t>(c, tag, input,
                                             make_state<eos_re_extra_t>(rho, T, 1.0e17, xn, junk), &p, &t);
                    cnt.other += 2;
                }
            }
        }
    }
}

} // namespace loom_eos_val

inline long validate_eos()
{
    using namespace loom_eos_val;

    LoomCheck c;
    c.module = "eos";
    std::mt19937_64 rng(20261018ULL);

    int n_test_react = 0;
    auto comps = make_compositions(rng, 30, &n_test_react);
    Counts cnt;

    // stock settings the parameter variants below override
    const int do_coulomb0 = helmholtz::do_coulomb;
    const int input_is_constant0 = helmholtz::input_is_constant;
    const Real rho_c0 = eos_rp::prad_limiter_rho_c;
    const Real delta_rho0 = eos_rp::prad_limiter_delta_rho;

    run_polys(c, rng);
    {
        mpl::EosTables t{};
        fill_eos_tables(&t);
        run_table_index(c, &t);
    }

    // the inputs file's settings: dense sweep plus edge cases
    sweep(c, "default", true, 1, comps, n_test_react, rng, cnt);
    edges(c, "default", comps, rng, cnt);

    // runtime-parameter variants (coarse sweeps, every 3rd composition)
    helmholtz::do_coulomb = 0;
    sweep(c, "no_coulomb", false, 3, comps, n_test_react, rng, cnt);
    helmholtz::do_coulomb = do_coulomb0;

    helmholtz::input_is_constant = 0;
    sweep(c, "input_not_constant", false, 3, comps, n_test_react, rng, cnt);
    edges(c, "input_not_constant", comps, rng, cnt);
    helmholtz::input_is_constant = input_is_constant0;

    // radiation pressure limiter (eos.prad_limiter_*; default off)
    for (Real rc : {1.0e2, 1.0e6, 1.0e8}) {
        eos_rp::prad_limiter_rho_c = rc;
        eos_rp::prad_limiter_delta_rho = 0.3 * rc;
        sweep(c, "prad_limiter", false, 3, comps, n_test_react, rng, cnt);
    }
    eos_rp::prad_limiter_rho_c = rho_c0;
    eos_rp::prad_limiter_delta_rho = delta_rho0;

    amrex::Print() << "LOOM_VALIDATE eos samples: compositions=" << comps.size()
                   << " (test_react " << n_test_react << ") rt_calls=" << cnt.rt
                   << " re_calls=" << cnt.re << " other_mode_calls=" << cnt.other
                   << " reset_inputs=" << cnt.reset << " pieces=" << cnt.pieces
                   << " composition_derivatives=" << cnt.xderivs << "\n";
    return loom_report(c);
}

#endif
