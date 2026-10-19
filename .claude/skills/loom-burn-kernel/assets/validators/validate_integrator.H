// Native bit-exact validation of vode.H + burner.H (host only; included by
// loom_validate.cpp under -DLOOM_VALIDATE_INTEGRATOR). See README.md.
//
// Run with integrator.use_jacobian_caching=0 so the stock CPU path (compiled
// with ALLOW_JACOBIAN_CACHING) takes the same branch as the HIP build; the
// validator forces it to 0 while it runs and restores it.
//
// 1. Pieces: normalize_abundances_burn; the integrator rhs()/jac() wrappers
//    on dvode_t states (every integrator_rp switch they read); dvjac with
//    the numerical Jacobian (jacobian_type = 2) on mid-burn-like states;
//    numerical_jac (integration/utils/numerical_jacobian.H).
// 2. Whole burns, exactly as test_react's do_react prepares each zone
//    (main.cpp rho/T log grids, setup_composition/get_xn, max(X, 1e-10),
//    normalize_abundances_burn, e = 0, T_fixed = -1, dt = unit_test.tmax):
//    stock burner(burn_state, dt) vs mpl::burner(...), every burn_t field
//    compared bit for bit, on n_cell = 16 and 32 grids, under a set of
//    integrator_rp variants (retry, Jacobian type, scaling, clipping, ...),
//    plus hand-picked hard zones.
// 3. Sensitivity: a 1-ulp perturbation of the Loom-side atol_spec must
//    produce mismatches (reported, not counted).
#ifndef VALIDATE_INTEGRATOR_H
#define VALIDATE_INTEGRATOR_H

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <functional>
#include <random>
#include <string>
#include <vector>

#include <extern_parameters.H>
#include <eos.H>             // stock (the kernel dir is last on -I)
#include <burn_type.H>
#include <burner.H>
#include <integrator.H>
#include <vode_type.H>
#include <vode_dvjac.H>
#include <integrator_rhs_strang.H>
#include <numerical_jacobian.H>
#include <react_util.H>
#include <loom_validate.H>

#include "validate_eos.H"     // EOS static_asserts
#include "validate_network.H" // burn_t layout static_asserts
#include "burner.H"           // translated (this directory)
#include <loom_params.H>       // host fill_* functions (unit_test/test_react)

// fill_integrator_params / fill_burn_params: <loom_params.H>

namespace loomval_integ {

using mpl::Real;

static_assert(INT_NEQS == mpl::neqs);
static_assert(VODE_LMAX == mpl::VODE_LMAX);
static_assert(UROUND == mpl::UROUND);
static_assert(CCMXJ == mpl::CCMXJ && HMIN == mpl::HMIN);
static_assert(vode_increase_change_factor == mpl::vode_increase_change_factor);
static_assert(vode_decrease_change_factor == mpl::vode_decrease_change_factor);
static_assert(vode_final_state_species_failure_tolerance_factor ==
              mpl::vode_final_state_species_failure_tolerance_factor);
static_assert(IERR_SUCCESS == mpl::IERR_SUCCESS && IERR_DT_UNDERFLOW == mpl::IERR_DT_UNDERFLOW &&
              IERR_TOO_MANY_STEPS == mpl::IERR_TOO_MANY_STEPS &&
              IERR_TOO_MUCH_ACCURACY_REQUESTED == mpl::IERR_TOO_MUCH_ACCURACY_REQUESTED &&
              IERR_CORRECTOR_CONVERGENCE == mpl::IERR_CORRECTOR_CONVERGENCE);
#ifdef SDC
#error "integrator translation is STRANG only"
#endif
#ifdef NSE
#error "integrator translation assumes no NSE"
#endif
#ifdef NONAKA_PLOT
#error "integrator translation assumes no NONAKA_PLOT"
#endif

using dvode14 = dvode_t<INT_NEQS>;

// ---------------------------------------------------------------------------
// stock-global snapshot / restore
// ---------------------------------------------------------------------------

struct Globals {
    Real rtol_spec, rtol_enuc, atol_spec, atol_enuc;
    Real retry_rtol_spec, retry_rtol_enuc, retry_atol_spec, retry_atol_enuc;
    int jacobian;
    bool use_burn_retry, retry_swap_jacobian, scale_system, call_eos_in_rhs, integrate_energy;
    bool use_number_densities, subtract_internal_energy, do_species_clip, renormalize_abundances;
    bool do_corrector_validation, correct_jacobian_for_const_e, use_jacobian_caching;
    int ode_max_steps;
    Real ode_max_dt, SMALL_X_SAFE, MAX_TEMP, react_boost, species_failure_tolerance, X_reject_buffer;
    Real small_x;
};

inline Globals save_globals()
{
    namespace ir = integrator_rp;
    return Globals{ir::rtol_spec, ir::rtol_enuc, ir::atol_spec, ir::atol_enuc,
                   ir::retry_rtol_spec, ir::retry_rtol_enuc, ir::retry_atol_spec, ir::retry_atol_enuc,
                   ir::jacobian,
                   static_cast<bool>(ir::use_burn_retry), static_cast<bool>(ir::retry_swap_jacobian),
                   static_cast<bool>(ir::scale_system), static_cast<bool>(ir::call_eos_in_rhs),
                   static_cast<bool>(ir::integrate_energy), static_cast<bool>(ir::use_number_densities),
                   static_cast<bool>(ir::subtract_internal_energy), static_cast<bool>(ir::do_species_clip),
                   static_cast<bool>(ir::renormalize_abundances), static_cast<bool>(ir::do_corrector_validation),
                   static_cast<bool>(ir::correct_jacobian_for_const_e), static_cast<bool>(ir::use_jacobian_caching),
                   ir::ode_max_steps, ir::ode_max_dt, ir::SMALL_X_SAFE, ir::MAX_TEMP, ir::react_boost,
                   ir::species_failure_tolerance, ir::X_reject_buffer, network_rp::small_x};
}

inline void restore_globals(const Globals& g)
{
    namespace ir = integrator_rp;
    ir::rtol_spec = g.rtol_spec; ir::rtol_enuc = g.rtol_enuc;
    ir::atol_spec = g.atol_spec; ir::atol_enuc = g.atol_enuc;
    ir::retry_rtol_spec = g.retry_rtol_spec; ir::retry_rtol_enuc = g.retry_rtol_enuc;
    ir::retry_atol_spec = g.retry_atol_spec; ir::retry_atol_enuc = g.retry_atol_enuc;
    ir::jacobian = g.jacobian;
    ir::use_burn_retry = g.use_burn_retry; ir::retry_swap_jacobian = g.retry_swap_jacobian;
    ir::scale_system = g.scale_system; ir::call_eos_in_rhs = g.call_eos_in_rhs;
    ir::integrate_energy = g.integrate_energy; ir::use_number_densities = g.use_number_densities;
    ir::subtract_internal_energy = g.subtract_internal_energy; ir::do_species_clip = g.do_species_clip;
    ir::renormalize_abundances = g.renormalize_abundances;
    ir::do_corrector_validation = g.do_corrector_validation;
    ir::correct_jacobian_for_const_e = g.correct_jacobian_for_const_e;
    ir::use_jacobian_caching = g.use_jacobian_caching;
    ir::ode_max_steps = g.ode_max_steps; ir::ode_max_dt = g.ode_max_dt;
    ir::SMALL_X_SAFE = g.SMALL_X_SAFE; ir::MAX_TEMP = g.MAX_TEMP; ir::react_boost = g.react_boost;
    ir::species_failure_tolerance = g.species_failure_tolerance;
    ir::X_reject_buffer = g.X_reject_buffer;
    network_rp::small_x = g.small_x;
}

// ---------------------------------------------------------------------------
// comparisons
// ---------------------------------------------------------------------------

inline void cb(LoomCheck& c, const std::string& w, double s, double m) { loom_check_bits(c, w.c_str(), s, m); }
inline void ci(LoomCheck& c, const std::string& w, long s, long m) { loom_check_int(c, w.c_str(), s, m); }

// every burn_t field
inline void compare_burn(LoomCheck& c, const std::string& tag, const burn_t& s, const mpl::burn_t& m)
{
    cb(c, tag + ".time", s.time, m.time);
    cb(c, tag + ".rho", s.rho, m.rho);
    cb(c, tag + ".T", s.T, m.T);
    cb(c, tag + ".e", s.e, m.e);
    for (int n = 0; n < NumSpec; ++n) { cb(c, tag + ".xn[" + std::to_string(n) + "]", s.xn[n], m.xn[n]); }
    cb(c, tag + ".e_scale", s.e_scale, m.e_scale);
    cb(c, tag + ".dedr", s.dedr, m.dedr);
    cb(c, tag + ".dedT", s.dedT, m.dedT);
    cb(c, tag + ".dedA", s.dedA, m.dedA);
    cb(c, tag + ".dedZ", s.dedZ, m.dedZ);
    ci(c, tag + ".i", s.i, m.i);
    ci(c, tag + ".j", s.j, m.j);
    ci(c, tag + ".k", s.k, m.k);
    for (int n = 0; n < SVAR; ++n) { cb(c, tag + ".y[" + std::to_string(n) + "]", s.y[n], m.y[n]); }
    cb(c, tag + ".rho_orig", s.rho_orig, m.rho_orig);
    cb(c, tag + ".rhoe_orig", s.rhoe_orig, m.rhoe_orig);
    for (int n = 0; n < SVAR; ++n) { cb(c, tag + ".ydot_a[" + std::to_string(n) + "]", s.ydot_a[n], m.ydot_a[n]); }
    ci(c, tag + ".sdc_iter", s.sdc_iter, m.sdc_iter);
    ci(c, tag + ".num_sdc_iters", s.num_sdc_iters, m.num_sdc_iters);
    cb(c, tag + ".T_fixed", s.T_fixed, m.T_fixed);
    cb(c, tag + ".cv", s.cv, m.cv);
    cb(c, tag + ".dx", s.dx, m.dx);
    cb(c, tag + ".mu", s.mu, m.mu);
    cb(c, tag + ".mu_e", s.mu_e, m.mu_e);
    cb(c, tag + ".y_e", s.y_e, m.y_e);
    cb(c, tag + ".eta", s.eta, m.eta);
    cb(c, tag + ".abar", s.abar, m.abar);
    cb(c, tag + ".zbar", s.zbar, m.zbar);
    ci(c, tag + ".n_rhs", s.n_rhs, m.n_rhs);
    ci(c, tag + ".n_jac", s.n_jac, m.n_jac);
    ci(c, tag + ".n_step", s.n_step, m.n_step);
    ci(c, tag + ".success", s.success ? 1 : 0, m.success ? 1 : 0);
    ci(c, tag + ".error_code", s.error_code, m.error_code);
}

// every dvode_t field (stock dvode_t<14> without jac_save/pivot)
inline void compare_dvode(LoomCheck& c, const std::string& tag, const dvode14& s, const mpl::dvode_t& m)
{
    cb(c, tag + ".CONP", s.CONP, m.CONP); cb(c, tag + ".CRATE", s.CRATE, m.CRATE);
    cb(c, tag + ".DRC", s.DRC, m.DRC); cb(c, tag + ".ETA", s.ETA, m.ETA);
    cb(c, tag + ".ETAMAX", s.ETAMAX, m.ETAMAX); cb(c, tag + ".H", s.H, m.H);
    cb(c, tag + ".HSCAL", s.HSCAL, m.HSCAL); cb(c, tag + ".PRL1", s.PRL1, m.PRL1);
    cb(c, tag + ".HMXI", s.HMXI, m.HMXI); cb(c, tag + ".RC", s.RC, m.RC);
    cb(c, tag + ".RL1", s.RL1, m.RL1); cb(c, tag + ".tn", s.tn, m.tn);
    ci(c, tag + ".n_rhs", s.n_rhs, m.n_rhs); ci(c, tag + ".n_jac", s.n_jac, m.n_jac);
    ci(c, tag + ".n_step", s.n_step, m.n_step);
    ci(c, tag + ".ICF", s.ICF, m.ICF); ci(c, tag + ".IPUP", s.IPUP, m.IPUP);
    ci(c, tag + ".JCUR", s.JCUR, m.JCUR); ci(c, tag + ".L", s.L, m.L);
    ci(c, tag + ".NEWH", s.NEWH, m.NEWH); ci(c, tag + ".NEWQ", s.NEWQ, m.NEWQ);
    ci(c, tag + ".NQ", s.NQ, m.NQ); ci(c, tag + ".NQWAIT", s.NQWAIT, m.NQWAIT);
    ci(c, tag + ".NSLJ", s.NSLJ, m.NSLJ); ci(c, tag + ".NSLP", s.NSLP, m.NSLP);
    ci(c, tag + ".jacobian_type", s.jacobian_type, m.jacobian_type);
    for (int i = 1; i <= VODE_LMAX; ++i) {
        cb(c, tag + ".el", s.el(i), m.el[i]);
        cb(c, tag + ".tau", s.tau(i), m.tau[i]);
    }
    for (int i = 1; i <= 5; ++i) { cb(c, tag + ".tq", s.tq(i), m.tq[i]); }
    cb(c, tag + ".rtol_spec", s.rtol_spec, m.rtol_spec); cb(c, tag + ".atol_spec", s.atol_spec, m.atol_spec);
    cb(c, tag + ".rtol_enuc", s.rtol_enuc, m.rtol_enuc); cb(c, tag + ".atol_enuc", s.atol_enuc, m.atol_enuc);
    cb(c, tag + ".t", s.t, m.t); cb(c, tag + ".tout", s.tout, m.tout);
    for (int i = 1; i <= INT_NEQS; ++i) {
        cb(c, tag + ".y", s.y(i), m.y[i]);
        cb(c, tag + ".ewt", s.ewt(i), m.ewt[i]);
        cb(c, tag + ".savf", s.savf(i), m.savf[i]);
        cb(c, tag + ".acor", s.acor(i), m.acor[i]);
        for (int j = 1; j <= INT_NEQS; ++j) { cb(c, tag + ".jac", s.jac(i, j), m.jac[MPL_JAC(i, j)]); }
        for (int j = 1; j <= VODE_LMAX; ++j) { cb(c, tag + ".yh", s.yh(i, j), m.yh[MPL_YH(i, j)]); }
    }
}

inline mpl::burn_t to_mpl(const burn_t& s)
{
    mpl::burn_t m;
    std::memcpy(&m, &s, sizeof m);  // identical layout (validate_network.H static_asserts)
    return m;
}

inline void copy_dvode(const dvode14& s, mpl::dvode_t* m)
{
    m->CONP = s.CONP; m->CRATE = s.CRATE; m->DRC = s.DRC; m->ETA = s.ETA; m->ETAMAX = s.ETAMAX;
    m->H = s.H; m->HSCAL = s.HSCAL; m->PRL1 = s.PRL1; m->HMXI = s.HMXI; m->RC = s.RC; m->RL1 = s.RL1;
    m->tn = s.tn; m->n_rhs = s.n_rhs; m->n_jac = s.n_jac; m->n_step = s.n_step;
    m->ICF = s.ICF; m->IPUP = s.IPUP; m->JCUR = s.JCUR; m->L = s.L; m->NEWH = s.NEWH; m->NEWQ = s.NEWQ;
    m->NQ = s.NQ; m->NQWAIT = s.NQWAIT; m->NSLJ = s.NSLJ; m->NSLP = s.NSLP; m->jacobian_type = s.jacobian_type;
    for (int i = 1; i <= VODE_LMAX; ++i) { m->el[i] = s.el(i); m->tau[i] = s.tau(i); }
    for (int i = 1; i <= 5; ++i) { m->tq[i] = s.tq(i); }
    m->rtol_spec = s.rtol_spec; m->atol_spec = s.atol_spec; m->rtol_enuc = s.rtol_enuc; m->atol_enuc = s.atol_enuc;
    m->t = s.t; m->tout = s.tout;
    for (int i = 1; i <= INT_NEQS; ++i) {
        m->y[i] = s.y(i); m->ewt[i] = s.ewt(i); m->savf[i] = s.savf(i); m->acor[i] = s.acor(i);
        for (int j = 1; j <= INT_NEQS; ++j) { m->jac[MPL_JAC(i, j)] = s.jac(i, j); }
        for (int j = 1; j <= VODE_LMAX; ++j) { m->yh[MPL_YH(i, j)] = s.yh(i, j); }
    }
}

struct Rng {
    std::mt19937_64 gen{0x766f6465ULL};
    Real uni(Real a, Real b) { return std::uniform_real_distribution<Real>(a, b)(gen); }
    int pick(int n) { return std::uniform_int_distribution<int>(0, n - 1)(gen); }
    Real logu(Real lo, Real hi) { return std::pow(10.0, uni(std::log10(lo), std::log10(hi))); }
};

// stock variants of the integrator_rp switches (each returns to the saved
// globals before the next one)
struct Variant {
    const char* name;
    std::function<void()> set;
};

inline std::vector<Variant> rhs_jac_variants()
{
    namespace ir = integrator_rp;
    return {
        {"default", [] {}},
        {"scale_system", [] { ir::scale_system = true; }},
        {"no_eos_in_rhs", [] { ir::call_eos_in_rhs = false; }},
        {"no_integrate_energy", [] { ir::integrate_energy = false; }},
        {"number_densities", [] { ir::use_number_densities = true; }},
        {"species_clip", [] { ir::do_species_clip = true; }},
        {"renormalize", [] { ir::renormalize_abundances = true; }},
        {"react_boost", [] { ir::react_boost = 3.5; }},
        {"no_const_e_corr", [] { ir::correct_jacobian_for_const_e = false; }},
        {"max_temp_low", [] { ir::MAX_TEMP = 2.0e9; }},
        {"all_on", [] { ir::scale_system = true; ir::do_species_clip = true; ir::renormalize_abundances = true;
                        ir::react_boost = 0.5; ir::SMALL_X_SAFE = 1.0e-12; }},
    };
}

// ---------------------------------------------------------------------------
// 1. pieces
// ---------------------------------------------------------------------------

inline void random_comp(Rng& rng, Real* x, bool allow_bad)
{
    Real sum = 0.0;
    for (int n = 0; n < NumSpec; ++n) {
        x[n] = (rng.pick(3) == 0) ? 0.0 : rng.logu(1.0e-14, 1.0);
        sum += x[n];
    }
    for (int n = 0; n < NumSpec; ++n) { x[n] /= sum; }
    if (allow_bad) {
        x[rng.pick(NumSpec)] = rng.uni(-0.05, 1.2);
        x[rng.pick(NumSpec)] = -rng.logu(1.0e-40, 1.0e-3);
    }
}

inline burn_t thermo_state(Rng& rng, Real rho, Real T, bool bad_comp)
{
    burn_t s;
    random_comp(rng, s.xn, false);
    s.rho = rho;
    s.T = T;
    s.T_fixed = -1.0;
    eos(eos_input_rt, s);
    if (bad_comp) {
        Real x[NumSpec];
        random_comp(rng, x, true);
        for (int n = 0; n < NumSpec; ++n) { s.xn[n] = x[n]; }
    }
    s.e_scale = (rng.pick(2) == 0) ? s.e : rng.logu(1.0e10, 1.0e20);
    return s;
}

inline void fill_dvode_random(Rng& rng, const burn_t& s, dvode14& v)
{
    v = dvode14{};
    for (int i = 1; i <= NumSpec; ++i) { v.y(i) = s.xn[i-1] * (1.0 + rng.uni(-1.0e-3, 1.0e-3)); }
    v.y(net_ienuc) = s.e / (integrator_rp::scale_system ? s.e_scale : 1.0) * (1.0 + rng.uni(-1.0e-6, 1.0e-6));
    v.t = 0.0;
    v.tout = 1.0e-5;
    v.tn = rng.uni(0.0, 1.0e-5);
    v.H = rng.logu(1.0e-14, 1.0e-6);
    v.RL1 = rng.uni(0.3, 1.0);
    v.n_step = rng.pick(100);
    v.n_rhs = rng.pick(1000);
    v.n_jac = rng.pick(50);
    v.jacobian_type = 2;
    v.rtol_spec = integrator_rp::rtol_spec; v.atol_spec = integrator_rp::atol_spec;
    v.rtol_enuc = integrator_rp::rtol_enuc; v.atol_enuc = integrator_rp::atol_enuc;
    for (int i = 1; i <= INT_NEQS; ++i) {
        v.ewt(i) = 1.0 / (v.rtol_spec * std::abs(v.y(i)) + (i <= NumSpec ? v.atol_spec : v.atol_enuc));
        v.acor(i) = rng.uni(-1.0, 1.0);
        for (int j = 1; j <= VODE_LMAX; ++j) { v.yh(i, j) = rng.uni(-1.0, 1.0); }
        for (int j = 1; j <= INT_NEQS; ++j) { v.jac(i, j) = rng.uni(-1.0, 1.0); }
    }
}

inline long g_rhs_zero_branch = 0;
inline long g_jac_zero_branch = 0;
inline long g_numjac_zero_branch = 0;

inline void run_pieces(LoomCheck& c)
{
    Rng rng;
    mpl::EosTables tab{};
    fill_eos_tables(&tab);
    const Globals g0 = save_globals();

    // normalize_abundances_burn
    for (int it = 0; it < 2000; ++it) {
        burn_t s;
        random_comp(rng, s.xn, it % 2 == 1);
        if (it % 7 == 0) { network_rp::small_x = rng.logu(1.0e-30, 1.0e-3); }
        mpl::burn_t m = to_mpl(s);
        mpl::IntegratorParams ip{};
        fill_integrator_params(&ip);
        normalize_abundances_burn(s);
        mpl::normalize_abundances_burn(&m, &ip);
        compare_burn(c, "normalize", s, m);
        restore_globals(g0);
    }

    const std::vector<Real> temps = {1.0e4, 1.0e6, 5.0e7, 3.0e8, 1.0e9, 3.0e9, 6.0e9, 1.0e10, 2.0e11};
    const std::vector<Real> rhos = {1.0e-2, 1.0e4, 1.0e6, 1.0e8, 1.0e10};

    for (const Variant& var : rhs_jac_variants()) {
        var.set();
        mpl::BurnParams bp{};
        fill_burn_params(&bp);
        for (Real T : temps) {
            for (Real rho : rhos) {
                for (int rep = 0; rep < 3; ++rep) {
                    const std::string tag = std::string("rhsjac.") + var.name;
                    burn_t s0 = thermo_state(rng, rho, T, rep == 2);
                    if (rep == 1) { s0.T_fixed = rng.uni(1.0e8, 5.0e9); }
                    dvode14 v0;
                    fill_dvode_random(rng, s0, v0);
                    const Real time = rng.uni(0.0, 1.0e-5);

                    // integrator rhs
                    {
                        burn_t s = s0; dvode14 v = v0;
                        mpl::burn_t m = to_mpl(s0); mpl::dvode_t mv{}; copy_dvode(v0, &mv);
                        RArray1D ydot;
                        for (int i = 1; i <= INT_NEQS; ++i) { ydot(i) = -7.0; }
                        Real mydot[mpl::neqs + 1];
                        for (int i = 0; i <= mpl::neqs; ++i) { mydot[i] = -7.0; }
                        ::rhs(time, s, v, ydot);
                        mpl::rhs(time, &m, &mv, mydot, &bp.integ, &bp.eos, &bp.net, tab);
                        if (s.T <= EOSData::mintemp || s.T >= integrator_rp::MAX_TEMP) { ++g_rhs_zero_branch; }
                        for (int i = 1; i <= INT_NEQS; ++i) { cb(c, tag + ".rhs.ydot", ydot(i), mydot[i]); }
                        compare_burn(c, tag + ".rhs.state", s, m);
                        compare_dvode(c, tag + ".rhs.vstate", v, mv);
                    }
                    // integrator jac (state as left by an rhs call, as in dvnlsd)
                    {
                        burn_t s = s0; dvode14 v = v0;
                        RArray1D ydot;
                        ::rhs(time, s, v, ydot);
                        mpl::burn_t m = to_mpl(s); mpl::dvode_t mv{}; copy_dvode(v, &mv);
                        ::jac(time, s, v, v.jac);
                        mpl::jac(time, &m, &mv, mv.jac, &bp.integ, &bp.eos, &bp.net, tab);
                        if (s.T <= EOSData::mintemp || s.T >= integrator_rp::MAX_TEMP) { ++g_jac_zero_branch; }
                        compare_burn(c, tag + ".jac.state", s, m);
                        compare_dvode(c, tag + ".jac.vstate", v, mv);
                    }
                    // dvjac: numerical (2) and analytic (1) Jacobian + LU
                    for (int jt = 1; jt <= 2; ++jt) {
                        burn_t s = s0; dvode14 v = v0;
                        v.jacobian_type = static_cast<short>(jt);
                        ::rhs(v.tn, s, v, v.savf);
                        mpl::burn_t m = to_mpl(s); mpl::dvode_t mv{}; copy_dvode(v, &mv);
                        int ier_s = -5;
                        int ier_m = -5;
                        ::dvjac(ier_s, s, v);
                        mpl::dvjac(ier_m, &m, &mv, &bp.integ, &bp.eos, &bp.net, tab);
                        ci(c, tag + ".dvjac.IERPJ", ier_s, ier_m);
                        compare_burn(c, tag + ".dvjac.state", s, m);
                        compare_dvode(c, tag + ".dvjac.vstate", v, mv);
                    }
                    // numerical_jac (integration/utils/numerical_jacobian.H)
                    {
                        burn_t s = s0;
                        if (rep == 2) { s.T_fixed = -1.0; }
                        jac_info_t ji{rng.logu(1.0e-14, 1.0e-4)};
                        JacNetArray2D J;
                        for (int i = 1; i <= INT_NEQS; ++i) {
                            for (int j = 1; j <= INT_NEQS; ++j) { J(i, j) = -3.0; }
                        }
                        Real mJ[MPL_JAC_SIZE];
                        for (int i = 0; i < MPL_JAC_SIZE; ++i) { mJ[i] = -3.0; }
                        mpl::burn_t m = to_mpl(s);
                        ::numerical_jac(s, ji, J);
                        mpl::numerical_jac(&m, ji.h, mJ, &bp.integ, &bp.eos, &bp.net, tab);
                        Real Tp = s.T + std::max(std::sqrt(U) * std::abs(s.T), 1.0);
                        if (Tp >= integrator_rp::MAX_TEMP) { ++g_numjac_zero_branch; }
                        for (int i = 1; i <= INT_NEQS; ++i) {
                            for (int j = 1; j <= INT_NEQS; ++j) {
                                cb(c, tag + ".numerical_jac", J(i, j), mJ[MPL_JAC(i, j)]);
                            }
                        }
                        compare_burn(c, tag + ".numerical_jac.state", s, m);
                    }
                }
            }
        }
        restore_globals(g0);
    }
}

// ---------------------------------------------------------------------------
// 2. whole burns
// ---------------------------------------------------------------------------

struct Zone {
    burn_t s;
    Real dt;
};

// test_react main.cpp initialization + react_zones.H do_react preparation
inline std::vector<Zone> test_react_zones(int n_cell)
{
    std::vector<Zone> zones;
    Real dlogrho = 0.0;
    Real dlogT = 0.0;
    if (n_cell > 1) {
        dlogrho = (std::log10(unit_test_rp::dens_max) - std::log10(unit_test_rp::dens_min)) / static_cast<amrex::Real>(n_cell - 1);
        dlogT = (std::log10(unit_test_rp::temp_max) - std::log10(unit_test_rp::temp_min)) / static_cast<amrex::Real>(n_cell - 1);
    }
    init_t comp_data = setup_composition(n_cell);
    for (int k = 0; k < n_cell; ++k) {
        for (int j = 0; j < n_cell; ++j) {
            for (int i = 0; i < n_cell; ++i) {
                const Real temp = std::pow(10.0_rt, (std::log10(unit_test_rp::temp_min) + static_cast<amrex::Real>(j)*dlogT));
                const Real dens = std::pow(10.0_rt, (std::log10(unit_test_rp::dens_min) + static_cast<amrex::Real>(i)*dlogrho));
                Real xn[NumSpec];
                get_xn(k, comp_data, xn, unit_test_rp::uniform_xn);

                burn_t burn_state;
                burn_state.rho = dens;
                burn_state.T = temp;
                for (int n = 0; n < NumSpec; ++n) {
                    burn_state.xn[n] = amrex::max(xn[n], 1.e-10_rt);
                }
                // normalize_abundances_burn(burn_state) is done (and checked)
                // in run_burn on both sides
                burn_state.e = 0.0_rt;
                burn_state.i = i;
                burn_state.j = j;
                burn_state.k = k;
                burn_state.T_fixed = -1.0_rt;
                zones.push_back(Zone{burn_state, unit_test_rp::tmax});
            }
        }
    }
    return zones;
}

inline std::vector<Zone> hard_zones()
{
    std::vector<Zone> z;
    auto add = [&z](Real rho, Real T, std::vector<std::pair<int, Real>> comp, Real dt, Real T_fixed) {
        burn_t s;
        for (int n = 0; n < NumSpec; ++n) { s.xn[n] = 1.0e-10; }
        for (auto& p : comp) { s.xn[p.first] = p.second; }
        s.rho = rho;
        s.T = T;
        s.e = 0.0;
        s.i = 99; s.j = 98; s.k = static_cast<int>(z.size());
        s.T_fixed = T_fixed;
        z.push_back(Zone{s, dt});
    };
    const Real tmax = unit_test_rp::tmax;
    // Network-specific hard zones (detonations, burning stages of specific
    // species, by this network's species indices) belong here; the template
    // keeps only edge cases that need nothing but species 0.
    add(1.0e5, 1.0e4, {{0, 1.0}}, tmax, -1.0);                    // T at the EOS floor
    add(1.0e8, 2.0e11, {{0, 1.0}}, tmax, -1.0);                   // T above MAX_TEMP
    add(1.0e-3, 1.0e9, {{0, 1.0}}, tmax, -1.0);                   // tiny density
    add(1.0e8, 3.0e9, {{0, 1.0}}, 0.0, -1.0);                     // dt = 0 (tout == t)
    add(1.0e8, 3.0e9, {{0, 1.0}}, -tmax, -1.0);                   // negative dt
    add(1.0e8, 3.0e9, {{0, 1.0}}, 1.0e-300, -1.0);                // tiny dt
    add(1.0e8, 3.0e9, {{0, 1.0}}, tmax, 4.0e9);                   // T_fixed > 0
    return z;
}

struct BurnStats {
    long zones{0};
    long stock_fail{0};
    long retried{0};
    std::vector<long> codes = std::vector<long>(12, 0);
};

inline void run_burn(LoomCheck& c, const std::string& tag, const Zone& z, const mpl::BurnParams* bp,
                     const mpl::EosTables& tab, BurnStats& st)
{
    burn_t s = z.s;
    mpl::burn_t m = to_mpl(z.s);

    normalize_abundances_burn(s);
    mpl::IntegratorParams ip = bp->integ;
    mpl::normalize_abundances_burn(&m, &ip);
    compare_burn(c, tag + ".normalized", s, m);

    // did the first attempt fail? (retry path taken)
    if (integrator_rp::use_burn_retry) {
        burn_t probe = s;
        actual_integrator(probe, z.dt);
        if (!probe.success) { ++st.retried; }
    }

    burner(s, z.dt);
    mpl::burner(&m, z.dt, bp, tab);
    compare_burn(c, tag, s, m);

    ++st.zones;
    if (!s.success) { ++st.stock_fail; }
    int code = s.error_code;
    st.codes[(code >= -8 && code <= 1) ? code + 8 : 10] += 1;
}

inline std::string stats_str(const BurnStats& st)
{
    std::string r = "zones=" + std::to_string(st.zones) + " stock_failures=" + std::to_string(st.stock_fail) +
                    " retried=" + std::to_string(st.retried) + " error_codes{";
    for (int i = 0; i <= 9; ++i) {
        if (st.codes[i] > 0) { r += " " + std::to_string(i - 8) + ":" + std::to_string(st.codes[i]); }
    }
    return r + " }";
}

inline std::vector<Variant> burn_variants()
{
    namespace ir = integrator_rp;
    return {
        {"jacobian2", [] { ir::jacobian = 2; }},
        {"retry_steps20", [] { ir::use_burn_retry = true; ir::ode_max_steps = 20; }},
        {"retry_steps150_noswap", [] { ir::use_burn_retry = true; ir::ode_max_steps = 150; ir::retry_swap_jacobian = false; }},
        {"retry_steps150_tols", [] { ir::use_burn_retry = true; ir::ode_max_steps = 150;
                                     ir::retry_rtol_spec = 1.0e-8; ir::retry_atol_spec = 1.0e-10;
                                     ir::retry_rtol_enuc = 1.0e-5; ir::retry_atol_enuc = 1.0e-4; }},
        {"retry_jac2_swap", [] { ir::use_burn_retry = true; ir::jacobian = 2; ir::ode_max_steps = 60; }},
        {"max_steps30", [] { ir::ode_max_steps = 30; }},
        {"scale_system", [] { ir::scale_system = true; }},
        {"no_subtract_e", [] { ir::subtract_internal_energy = false; }},
        {"no_eos_in_rhs", [] { ir::call_eos_in_rhs = false; }},
        {"no_integrate_energy", [] { ir::integrate_energy = false; }},
        {"react_boost", [] { ir::react_boost = 10.0; }},
        {"clip_renorm", [] { ir::do_species_clip = true; ir::renormalize_abundances = true; }},
        {"no_corrector_validation", [] { ir::do_corrector_validation = false; }},
        {"no_const_e_corr", [] { ir::correct_jacobian_for_const_e = false; }},
        {"number_densities", [] { ir::use_number_densities = true; }},
        // dt cap active but bounded: >= 256 steps per zone (ode_max_dt = 1e-8
        // with tmax = 1e-3 meant >= 1e5 steps per zone and ~10 min of silence
        // under Rosetta; max_steps30 covers the step-limit failure path)
        {"ode_max_dt", [] { ir::ode_max_dt = unit_test_rp::tmax / 256.0; ir::ode_max_steps = 4000; }},
        {"max_temp_3e9", [] { ir::MAX_TEMP = 3.0e9; }},
        {"x_reject_buffer", [] { ir::X_reject_buffer = 1.0e4; ir::species_failure_tolerance = 1.0e-4; }},
        {"loose_tols", [] { ir::rtol_spec = 1.0e-6; ir::atol_spec = 1.0e-6; ir::rtol_enuc = 1.0e-4; }},
        {"tight_tols", [] { ir::rtol_spec = 1.0e-15; ir::atol_spec = 1.0e-15; ir::rtol_enuc = 1.0e-12;
                            ir::atol_enuc = 1.0e-12; ir::ode_max_steps = 3000; }},
        {"small_x", [] { network_rp::small_x = 1.0e-8; }},
        {"too_much_accuracy", [] { ir::rtol_spec = 1.0e-20; ir::atol_spec = 1.0e-30; }},
        {"too_much_accuracy_late", [] { ir::rtol_spec = 1.0e-17; ir::atol_spec = 1.0e-22;
                                        ir::rtol_enuc = 1.0e-17; ir::atol_enuc = 1.0e-22; ir::ode_max_steps = 2000; }},
        {"boost_1e6_nocv", [] { ir::react_boost = 1.0e6; ir::do_corrector_validation = false; ir::ode_max_steps = 2000; }},
        {"boost_1e6_retry", [] { ir::react_boost = 1.0e6; ir::use_burn_retry = true; ir::ode_max_steps = 500;
                                 ir::retry_rtol_spec = 1.0e-6; ir::retry_atol_enuc = 1.0e-2; }},
    };
}

inline long validate_integrator()
{
    LoomCheck c{"integrator"};
    LoomCheck cp{"integrator_pieces"};

    const Globals g0 = save_globals();
    integrator_rp::use_jacobian_caching = false;   // the HIP (no caching) path
    const Globals gbase = save_globals();

    mpl::EosTables tab{};
    fill_eos_tables(&tab);

    // 1. pieces
    run_pieces(cp);
    amrex::Print() << "LOOM_VALIDATE integrator_pieces rhs_zero_branch=" << g_rhs_zero_branch
                   << " jac_zero_branch=" << g_jac_zero_branch
                   << " numjac_T_out_of_range=" << g_numjac_zero_branch << "\n";
    restore_globals(gbase);

    // 2. whole burns, default parameters
    const std::vector<Zone> z16 = test_react_zones(16);
    const std::vector<Zone> z32 = test_react_zones(32);
    const std::vector<Zone> zh = hard_zones();
    {
        mpl::BurnParams bp{};
        fill_burn_params(&bp);
        BurnStats st16, st32, sth;
        for (const Zone& z : z16) { run_burn(c, "burn16", z, &bp, tab, st16); }
        amrex::Print() << "LOOM_VALIDATE integrator burn n_cell=16 default " << stats_str(st16) << "\n";
        for (const Zone& z : z32) { run_burn(c, "burn32", z, &bp, tab, st32); }
        amrex::Print() << "LOOM_VALIDATE integrator burn n_cell=32 default " << stats_str(st32) << "\n";
        for (const Zone& z : zh) { run_burn(c, "burn_hard", z, &bp, tab, sth); }
        amrex::Print() << "LOOM_VALIDATE integrator burn hard_zones default " << stats_str(sth) << "\n";
    }

    // variants on the n_cell = 16 grid (every 3rd zone) + all hard zones
    for (const Variant& var : burn_variants()) {
        // progress line first: some variants take minutes, which otherwise
        // looks like a hang (E2E runs 2026-10-19)
        amrex::Print() << "LOOM_VALIDATE integrator burn variant=" << var.name << " start\n";
        restore_globals(gbase);
        var.set();
        mpl::BurnParams bp{};
        fill_burn_params(&bp);
        BurnStats st;
        const std::string tag = std::string("burn.") + var.name;
        for (std::size_t iz = 0; iz < z16.size(); iz += 3) { run_burn(c, tag, z16[iz], &bp, tab, st); }
        for (const Zone& z : zh) { run_burn(c, tag + ".hard", z, &bp, tab, st); }
        amrex::Print() << "LOOM_VALIDATE integrator burn variant=" << var.name << " " << stats_str(st) << "\n";
    }
    restore_globals(gbase);

    // 3. sensitivity: perturb one Loom-side parameter by 1 ulp
    {
        mpl::BurnParams bp{};
        fill_burn_params(&bp);
        bp.integ.integrator_atol_spec = std::nextafter(bp.integ.integrator_atol_spec, 1.0);
        LoomCheck cs{"integrator_sensitivity"};
        cs.reported = 20;  // expected mismatches: do not print them
        BurnStats st;
        for (std::size_t iz = 0; iz < z16.size(); iz += 37) { run_burn(cs, "sens", z16[iz], &bp, tab, st); }
        amrex::Print() << "LOOM_VALIDATE integrator_sensitivity (atol_spec +1ulp, expected nonzero) checks="
                       << cs.checks << " mismatches=" << cs.mismatches << "\n";
        if (cs.mismatches == 0) {
            amrex::Print() << "LOOM_MISMATCH integrator_sensitivity: perturbation went undetected\n";
            c.mismatches += 1;
        }
    }

    restore_globals(g0);
    return loom_report(cp) + loom_report(c);
}

} // namespace loomval_integ

using loomval_integ::validate_integrator;

#endif
