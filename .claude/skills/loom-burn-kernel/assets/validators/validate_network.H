// Native bit-exact validation of the translated network.H + network_rhs.H
// against the stock network (host only; included by loom_validate.cpp under
// -DLOOM_VALIDATE_NETWORK, and by validate_integrator.H for the burn_t layout
// static_asserts). Network-neutral template; per-port hooks live in
// validate_port.H.
//
// stock RHS::rhs<neqs> / RHS::jac vs mpl::actual_rhs / mpl::actual_jac on
// every element, and RHS::dgefa / dgesl vs rhs_dgefa / rhs_dgesl on
// I - gamma*J, over states built with the stock EOS (T 1e4..1e16,
// rho 1e2..1e12, compositions built from NumSpec, random states), for every
// parameter variant from validate_port.H and enable_debye_huckel_skip = 0/1.
//
// Snapshot of the triple_alpha_plus_cago port's validator (2026-10-18), made
// generic. A port may add checks of its own (aprox13 also compared
// construct_rate and the rate tables), but should keep these.
#ifndef VALIDATE_NETWORK_H
#define VALIDATE_NETWORK_H

#include <cmath>
#include <cstddef>
#include <cstring>
#include <random>
#include <vector>

#include <extern_parameters.H>
#include <eos.H>
#include <burn_type.H>
#include <rhs.H>
#include <integrator_data.H>
#include <loom_validate.H>

#include "network.H"

// ---------------------------------------------------------------------------
// burn_t layout: stock vs mpl (validators memcpy between the two). These are
// the fields of the configurations gen_network_data.py supports (no aux
// species, no SDC); PORT: add an LV_OFF line for any field burn.H adds.
// ---------------------------------------------------------------------------
static_assert(sizeof(burn_t) == sizeof(mpl::burn_t));
#define LV_OFF(f) static_assert(offsetof(burn_t, f) == offsetof(mpl::burn_t, f), "burn_t." #f);
LV_OFF(time) LV_OFF(rho) LV_OFF(T) LV_OFF(e) LV_OFF(xn) LV_OFF(e_scale) LV_OFF(dedr)
LV_OFF(dedT) LV_OFF(dedA) LV_OFF(dedZ) LV_OFF(i) LV_OFF(j) LV_OFF(k) LV_OFF(y)
LV_OFF(rho_orig) LV_OFF(rhoe_orig) LV_OFF(ydot_a) LV_OFF(sdc_iter) LV_OFF(num_sdc_iters)
LV_OFF(T_fixed) LV_OFF(cv) LV_OFF(dx) LV_OFF(mu) LV_OFF(mu_e) LV_OFF(y_e) LV_OFF(eta)
LV_OFF(abar) LV_OFF(zbar) LV_OFF(n_rhs) LV_OFF(n_jac) LV_OFF(n_step) LV_OFF(success)
LV_OFF(error_code)
#undef LV_OFF
static_assert(NumSpec == mpl::NumSpec && neqs == mpl::neqs && net_ienuc == mpl::net_ienuc &&
              INT_NEQS == mpl::neqs);
static_assert(Rates::NumRates == mpl::Rates::NumRates);

#ifdef LOOM_VALIDATE_NETWORK
#include <loom_params.H>
#include "validate_port.H"

namespace loomval_net {

using mpl::Real;

inline burn_t make_state(Real rho, Real T, const std::vector<Real>& X)
{
    burn_t s;
    s.rho = rho;
    s.T = T;
    for (int n = 0; n < NumSpec; ++n) { s.xn[n] = X[n]; }
    eos(eos_input_rt, s);
    return s;
}

// noinline: see validate_rates.H (GCC constant-folding of libm calls)
__attribute__((noinline)) void one(LoomCheck& c, const burn_t& s, const mpl::NetParams* np)
{
    mpl::burn_t m;
    std::memcpy(&m, &s, sizeof m);

    amrex::Array1D<Real, 1, neqs> ydot;
    for (int i = 1; i <= neqs; ++i) { ydot(i) = 0.0; }
    RHS::rhs<neqs>(s, ydot);
    mpl::Array1D<Real, 1, mpl::neqs> ly{};
    mpl::actual_rhs(m, *np, ly);
    for (int i = 1; i <= neqs; ++i) {
        loom_check_bits(c, "rhs", ydot(i), ly(i));
    }

    JacNetArray2D jac;
    jac.zero();
    RHS::jac(s, jac);
    mpl::JacNetArray2D lj{};
    mpl::actual_jac(m, *np, lj);
    for (int i = 1; i <= neqs; ++i) {
        for (int j = 1; j <= neqs; ++j) {
            loom_check_bits(c, "jac", jac(i, j), lj(i, j));
        }
    }

    // LU on I - gamma*J (the shape VODE factors), then a solve
    for (Real gamma : {1.0e-12, 1.0e-6, 1.0e-2}) {
        RArray2D a;
        mpl::JacNetArray2D la{};
        for (int i = 1; i <= neqs; ++i) {
            for (int j = 1; j <= neqs; ++j) {
                a(i, j) = (i == j ? 1.0 : 0.0) - gamma * jac(i, j);
                la(i, j) = a(i, j);
            }
        }
        const int info_s = RHS::dgefa(a);
        const int info_l = mpl::rhs_dgefa(la);
        loom_check_int(c, "dgefa.info", info_s, info_l);
        for (int i = 1; i <= neqs; ++i) {
            for (int j = 1; j <= neqs; ++j) {
                loom_check_bits(c, "dgefa", a(i, j), la(i, j));
            }
        }
        RArray1D b;
        mpl::Array1D<Real, 1, mpl::neqs> lb{};
        for (int i = 1; i <= neqs; ++i) { b(i) = ydot(i); lb(i) = ydot(i); }
        RHS::dgesl(a, b);
        mpl::rhs_dgesl(la, lb);
        for (int i = 1; i <= neqs; ++i) { loom_check_bits(c, "dgesl", b(i), lb(i)); }
    }
}

// one-hot per species, all pairs, uniform, near-pure, random (12 decades)
inline std::vector<std::vector<Real>> compositions(std::mt19937_64& rng, int n_random)
{
    std::vector<std::vector<Real>> comps;
    for (int a = 0; a < NumSpec; ++a) {
        std::vector<Real> x(NumSpec, 0.0);
        x[a] = 1.0;
        comps.push_back(x);
        for (int b = a + 1; b < NumSpec; ++b) {
            std::vector<Real> y(NumSpec, 0.0);
            y[a] = 0.7;
            y[b] = 0.3;
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
        for (auto& v : x) { v = std::pow(10.0, -12.0 * u01(rng)); s += v; }
        for (auto& v : x) { v /= s; }
        comps.push_back(x);
    }
    return comps;
}

} // namespace loomval_net

inline long validate_network()
{
    using namespace loomval_net;
    LoomCheck c;
    c.module = "network";

    std::mt19937_64 rng(4242);
    std::uniform_real_distribution<Real> u01(0.0, 1.0);
    const auto comps = compositions(rng, 10);
    // keep the grid affordable for large networks (aprox13: 105 pair compositions)
    const Real dlr = comps.size() > 40 ? 1.0 : 0.5;
    const Real dlt = comps.size() > 40 ? 0.4 : 0.2;

    const bool save_skip = screening_rp::enable_debye_huckel_skip;
    loomval_port::save_params();
    long nstates = 0;
    for (int v = 0; v < loomval_port::num_param_variants; ++v) {
        for (int skip = 0; skip <= 1; ++skip) {
            loomval_port::apply_param_variant(v);
            screening_rp::enable_debye_huckel_skip = skip;
            mpl::NetParams np{};
            fill_net_params(&np);
            for (const auto& X : comps) {
                for (Real lr = 2.0; lr <= 12.0; lr += dlr) {
                    for (Real lt = 4.0; lt <= 16.0; lt += dlt) {
                        one(c, make_state(std::pow(10.0, lr), std::pow(10.0, lt), X), &np);
                        ++nstates;
                    }
                }
            }
            for (int k = 0; k < 3000; ++k) {
                const auto& X = comps[k % comps.size()];
                one(c, make_state(std::pow(10.0, 2.0 + 10.0 * u01(rng)),
                                  std::pow(10.0, 4.0 + 12.0 * u01(rng)), X), &np);
                ++nstates;
            }
        }
    }
    loomval_port::restore_params();
    screening_rp::enable_debye_huckel_skip = save_skip;
    amrex::Print() << "LOOM_VALIDATE network states=" << nstates << "\n";
    return loom_report(c);
}
#endif // LOOM_VALIDATE_NETWORK

#endif
