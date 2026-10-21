// Screening for the Loom burn kernel (see README.md): screening/screen.H
// (screen5 configuration) and screening/screen_data.H, for both number types
// the burn uses:
//   *_r  number_t = amrex::Real      (networks/rhs.H rhs())
//   *_d  number_t = autodiff::dual   (networks/rhs.H jac(); temp_t = dual)
// The _r functions mirror the stock arithmetic exactly; the _d functions are
// the dual.H call sequences autodiff evaluates (dual.H header, section 3).
//
// Screen factors: networks/rhs.H apply_screening evaluates
// scrn::calculate_screen_factor (gcem::pow / gcem::log) as a `constexpr auto`,
// i.e. at compile time, in both the CPU and the HIP build. Nothing of it runs
// on the device, so screen_factor<Z1, A1, Z2, A2>() returns those
// compile-time values, generated from the stock constexpr evaluation into
// network_data.H (net_screen_*; round-trip literals) for every pair the
// network's apply_screening uses. gen_network_data.py and validate_rates.H
// check each field bit for bit against the stock constexpr evaluation.
// Likewise gamma_e_constants in fill_plasma_state (constexpr, gcem::pow).
#ifndef MPL_SCREEN_H
#define MPL_SCREEN_H

#include "support.H"
#include "dual.H"
#include "network_data.H"
#include "burn.H"
#include "rates.H"

#if MPL_GEN_SCREENING && !defined(MPL_GEN_SCREEN_METHOD_screen5)
#error "screen.H translates SCREEN_METHOD=screen5 only (network_data.H was generated for another method)"
#endif

namespace mpl {

// from screening/screen.H:54 plasma_state_t
template <typename T>
struct plasma_state_t {
    T qlam0z;
    T taufac;
    T aa;
    T temp;
    T zbar;
    T z2bar;
    T n_e;
    T gamma_e_fac;
};

// from screening/screen_data.H:13 screen_factors_t (SCREEN_METHOD_screen5
// fields; same order as the stock class)
struct screen_factors_t {
    Real z1;
    Real z2;
    Real a1;
    Real a2;
    Real zs53;
    Real z1_53;
    Real z2_53;
    Real zs13;
    Real zs13inv;
    Real zhat;
    Real zhat2;
    Real lzav;
    Real aznut;
};

// from screening/screen_data.H:80 calculate_screen_factor, as folded at
// compile time by networks/rhs.H apply_screening (`constexpr auto scn_fac =
// scrn::calculate_screen_factor(Z1, A1, Z2, A2)`), for the pairs this
// network screens: the values come from network_data.H (generated from the
// stock constexpr evaluation and checked bit for bit against it).
template <int Z1, int A1, int Z2, int A2>
MPL_DEVICE screen_factors_t screen_factor()
{
    constexpr int pair = net_screen_pair(Z1, A1, Z2, A2);
    static_assert(pair >= 0, "screen_factor: pair not screened by this network (network_data.H)");

    constexpr Real z1 = net_screen_z1(pair);
    constexpr Real z2 = net_screen_z2(pair);
    constexpr Real a1 = net_screen_a1(pair);
    constexpr Real a2 = net_screen_a2(pair);
    constexpr Real zs53 = net_screen_zs53(pair);
    constexpr Real z1_53 = net_screen_z1_53(pair);
    constexpr Real z2_53 = net_screen_z2_53(pair);
    constexpr Real zs13 = net_screen_zs13(pair);
    constexpr Real zs13inv = net_screen_zs13inv(pair);
    constexpr Real zhat = net_screen_zhat(pair);
    constexpr Real zhat2 = net_screen_zhat2(pair);
    constexpr Real lzav = net_screen_lzav(pair);
    constexpr Real aznut = net_screen_aznut(pair);

    screen_factors_t r{z1, z2, a1, a2, zs53, z1_53, z2_53, zs13, zs13inv,
                       zhat, zhat2, lzav, aznut};
    return r;
}

// from screening/screen.H:115 fill_plasma_state (number_t = temp_t = Real)
// y(1..NumSpec) = rhs_state.y
MPL_DEVICE void fill_plasma_state_r(plasma_state_t<Real>* state, Real temp,
                                    const Real dens, const NetYArray1D& y)
{
    Real sum = 0.0;
    for (int n = 1; n <= NumSpec; n++) {
        sum += y(n);
    }
    Real abar = 1.0 / sum;
    Real ytot = sum;

    Real m1 = 0.0;
    Real m2 = 0.0;
    for (int n = 1; n <= NumSpec; n++) {
        m1 += zion[n-1]*y(n);
        m2 += zion[n-1]*zion[n-1]*y(n);
    }

    Real temp_ad = temp;

    const Real co2 = (1.0/3.0) * 4.248719e3;

    Real zbar = m1 * abar;
    Real z2bar = m2 * abar;

    Real rr = dens * ytot;
    Real tempi = 1.0 / temp_ad;

    Real pp = math::sqrt(rr*tempi*(z2bar + zbar));

    state->qlam0z = 1.88e8 * tempi * pp;

    state->taufac = co2 * math::cbrt(tempi);

    Real xni = math::cbrt(rr * zbar);

    state->aa = 2.27493e5 * tempi * xni;

    state->temp = temp_ad;
    state->zbar = zbar;
    state->z2bar = z2bar;

    // static_cast<amrex::Real>(C::n_A)
    constexpr Real n_A = 6.02214076e+23;
    state->n_e = zbar * rr * n_A;

    // C::q_e * C::q_e * gcem::pow(4.0_rt * M_PI / 3.0_rt, 1.0_rt / 3.0_rt) / C::k_B
    constexpr Real gamma_e_constants = 0.0026936538176192486;
    state->gamma_e_fac = gamma_e_constants * math::cbrt(state->n_e);
}

// from screening/screen.H:115 fill_plasma_state (number_t = temp_t = autodiff::dual)
MPL_DEVICE void fill_plasma_state_d(plasma_state_t<Dual>* state, Dual temp,
                                    const Real dens, const NetYArray1D& y)
{
    Real sum = 0.0;
    for (int n = 1; n <= NumSpec; n++) {
        sum += y(n);
    }
    Real abar = 1.0 / sum;
    Real ytot = sum;

    Dual m1 = dual_const(0.0);
    Dual m2 = dual_const(0.0);
    for (int n = 1; n <= NumSpec; n++) {
        m1 = dual_add(m1, dual_const(zion[n-1]*y(n)));
        m2 = dual_add(m2, dual_const(zion[n-1]*zion[n-1]*y(n)));
    }

    Dual temp_ad = temp;

    const Real co2 = (1.0/3.0) * 4.248719e3;

    Dual zbar = dual_mul_num(m1, abar);
    Dual z2bar = dual_mul_num(m2, abar);

    Real rr = dens * ytot;
    Dual tempi = dual_div(dual_const(1.0), temp_ad);

    Dual pp = dual_sqrt(dual_mul(dual_mul_num(tempi, rr), dual_add(z2bar, zbar)));

    state->qlam0z = dual_mul(dual_mul_num(tempi, 1.88e8), pp);

    state->taufac = dual_mul(dual_const(co2), dual_cbrt(tempi));

    Dual xni = dual_cbrt(dual_mul_num(zbar, rr));

    state->aa = dual_mul(dual_mul_num(tempi, 2.27493e5), xni);

    state->temp = temp_ad;
    state->zbar = zbar;
    state->z2bar = z2bar;

    constexpr Real n_A = 6.02214076e+23;
    state->n_e = dual_mul_num(dual_mul_num(zbar, rr), n_A);

    constexpr Real gamma_e_constants = 0.0026936538176192486;
    state->gamma_e_fac = dual_mul(dual_const(gamma_e_constants), dual_cbrt(state->n_e));
}

// from screening/screen.H:187 debye_huckel (number_t = Real)
MPL_DEVICE Real debye_huckel_r(const plasma_state_t<Real>& state,
                               const screen_factors_t& scn_fac)
{
    Real z1z2 = scn_fac.z1 * scn_fac.z2;
    Real z_ratio = state.z2bar / state.zbar;

    Real Gamma_e = state.gamma_e_fac / state.temp;

    Real h_DH = z1z2 * math::sqrt(3 * math::powi<3>(Gamma_e) * z_ratio);

    return math::max(h_DH, 0.0);
}

// from screening/screen.H:187 debye_huckel (number_t = autodiff::dual)
MPL_DEVICE Dual debye_huckel_d(const plasma_state_t<Dual>& state,
                               const screen_factors_t& scn_fac)
{
    Real z1z2 = scn_fac.z1 * scn_fac.z2;
    Dual z_ratio = dual_div(state.z2bar, state.zbar);

    Dual Gamma_e = dual_div(state.gamma_e_fac, state.temp);

    // 3 * powi<3>(Gamma_e) is NDM(3, dual): the int 3 converts to 3.0
    Dual h_DH = dual_mul(dual_const(z1z2),
                         dual_sqrt(dual_mul(dual_mul_num(dual_powi<3>(Gamma_e), 3.0), z_ratio)));

    return dual_max(h_DH, dual_const(0.0));
}

// from screening/screen.H:213 actual_screen5 (number_t = Real)
MPL_DEVICE Real actual_screen5_r(const plasma_state_t<Real>& state,
                                 const screen_factors_t& scn_fac)
{
    constexpr Real fact    = 1.25992104989487e0;
    constexpr Real gamefx  = 0.3e0;
    constexpr Real gamefs  = 0.8e0;
    constexpr Real inv_dg = 1.0 / (gamefs - gamefx);

    const Real z1 = scn_fac.z1;
    const Real z2 = scn_fac.z2;

    const Real bb = z1 * z2;
    Real gamp = state.aa;

    Real qq = fact * bb * scn_fac.zs13inv;

    Real gamef = qq * gamp;

    Real tau12 = state.taufac * scn_fac.aznut;

    Real alph12 = gamef / tau12;

    if (alph12 > 1.6) {
        alph12   = 1.6e0;

        gamef    = 1.6e0 * tau12;

        qq = scn_fac.zs13/(fact * bb);
        gamp = gamef * qq;
    }

    const Real h12w = bb * state.qlam0z;

    Real h12 = h12w;

    if (gamef > gamefx) {

        const Real gamp14 = math::sqrt(math::sqrt(gamp));

        const Real cc = 0.896434e0 * gamp * scn_fac.zhat
            - 3.44740e0 * gamp14 * scn_fac.zhat2
            - 0.5551e0 * (math::log(gamp) + scn_fac.lzav)
            - 2.996e0;

        const Real a3 = alph12 * alph12 * alph12;

        const Real rr = (5.0/32.0) - alph12*(0.014e0 + 0.0128e0*alph12);

        const Real ss = tau12*rr;

        const Real tt = -0.0098e0 + 0.0048e0*alph12;

        const Real uu = 0.0055e0 + alph12*tt;

        const Real vv = gamef * alph12 * uu;

        h12 = cc - a3 * (ss + vv);

        const Real xlgfac = math::max(1.0 - 0.0562e0*a3, 0.77);

        h12 += math::log(xlgfac);

        if (gamef <= gamefs) {
            const Real weight1 =  inv_dg * (gamefs - gamef);
            const Real weight2 = inv_dg * (gamef - gamefx);

            h12 = h12w * weight1 + h12 * weight2;
        }
    }

    return math::max(h12, 0.0);
}

// from screening/screen.H:213 actual_screen5 (number_t = autodiff::dual)
MPL_DEVICE Dual actual_screen5_d(const plasma_state_t<Dual>& state,
                                 const screen_factors_t& scn_fac)
{
    constexpr Real fact    = 1.25992104989487e0;
    constexpr Real gamefx  = 0.3e0;
    constexpr Real gamefs  = 0.8e0;
    constexpr Real inv_dg = 1.0 / (gamefs - gamefx);

    const Real z1 = scn_fac.z1;
    const Real z2 = scn_fac.z2;

    const Real bb = z1 * z2;
    Dual gamp = state.aa;

    Real qq = fact * bb * scn_fac.zs13inv;

    Dual gamef = dual_mul_num(gamp, qq);

    Dual tau12 = dual_mul_num(state.taufac, scn_fac.aznut);

    Dual alph12 = dual_div(gamef, tau12);

    if (alph12.val > 1.6) {
        alph12   = dual_const(1.6e0);

        gamef    = dual_mul_num(tau12, 1.6e0);

        qq = scn_fac.zs13/(fact * bb);
        gamp = dual_mul_num(gamef, qq);
    }

    const Dual h12w = dual_mul_num(state.qlam0z, bb);

    Dual h12 = h12w;

    if (gamef.val > gamefx) {

        const Dual gamp14 = dual_sqrt(dual_sqrt(gamp));

        Dual s = dual_mul_num(dual_mul_num(gamp, 0.896434e0), scn_fac.zhat);
        s = dual_sub(s, dual_mul_num(dual_mul_num(gamp14, 3.44740e0), scn_fac.zhat2));
        s = dual_sub(s, dual_mul(dual_const(0.5551e0), dual_add_num(dual_log(gamp), scn_fac.lzav)));
        const Dual cc = dual_add_num(s, -2.996e0);

        const Dual a3 = dual_mul(dual_mul(alph12, alph12), alph12);

        const Dual rr = dual_sub(dual_const(5.0/32.0),
                                 dual_mul(alph12, dual_add_ndm(dual_const(0.014e0), 0.0128e0, alph12)));

        const Dual ss = dual_mul(tau12, rr);

        const Dual tt = dual_add_ndm(dual_const(-0.0098e0), 0.0048e0, alph12);

        const Dual uu = dual_add(dual_const(0.0055e0), dual_mul(alph12, tt));

        const Dual vv = dual_mul(dual_mul(gamef, alph12), uu);

        h12 = dual_sub(cc, dual_mul(a3, dual_add(ss, vv)));

        const Dual xlgfac = dual_max(dual_add_ndm(dual_const(1.0), -0.0562e0, a3), dual_const(0.77));

        h12 = dual_add(h12, dual_log(xlgfac));

        if (gamef.val <= gamefs) {
            const Dual weight1 = dual_mul(dual_const(inv_dg), dual_sub(dual_const(gamefs), gamef));
            const Dual weight2 = dual_mul(dual_const(inv_dg), dual_add_num(gamef, -gamefx));

            h12 = dual_add(dual_mul(h12w, weight1), dual_mul(h12, weight2));
        }
    }

    return dual_max(h12, dual_const(0.0));
}

// from screening/screen.H:656 actual_log_screen (number_t = Real)
MPL_DEVICE Real actual_log_screen_r(const plasma_state_t<Real>& state,
                                    const screen_factors_t& scn_fac,
                                    const RatesParams& params)
{
    Real log_scor = 0.0;
    bool done = false;
    if (params.screening_enable_debye_huckel_skip) {
        log_scor = debye_huckel_r(state, scn_fac);
        if (log_scor <= params.screening_debye_huckel_skip_threshold) {
            done = true;
        }
    }
    if (!done) {
        log_scor = actual_screen5_r(state, scn_fac);
    }
    return log_scor;
}

// from screening/screen.H:656 actual_log_screen (number_t = autodiff::dual)
MPL_DEVICE Dual actual_log_screen_d(const plasma_state_t<Dual>& state,
                                    const screen_factors_t& scn_fac,
                                    const RatesParams& params)
{
    Dual log_scor = dual_const(0.0);
    bool done = false;
    if (params.screening_enable_debye_huckel_skip) {
        log_scor = debye_huckel_d(state, scn_fac);
        if (log_scor.val <= params.screening_debye_huckel_skip_threshold) {
            done = true;
        }
    }
    if (!done) {
        log_scor = actual_screen5_d(state, scn_fac);
    }
    return log_scor;
}

// from screening/screen.H:708 actual_screen (number_t = Real)
MPL_DEVICE Real actual_screen_r(const plasma_state_t<Real>& state,
                                const screen_factors_t& scn_fac,
                                const RatesParams& params)
{
    Real log_scor = actual_log_screen_r(state, scn_fac, params);
    Real result = 1.0;
    if (!(log_scor <= 0.0)) {
        constexpr Real log_scor_max = 300.0;
        log_scor = math::min(log_scor, log_scor_max);

        result = math::exp(log_scor);
    }
    return result;
}

// from screening/screen.H:708 actual_screen (number_t = autodiff::dual)
MPL_DEVICE Dual actual_screen_d(const plasma_state_t<Dual>& state,
                                const screen_factors_t& scn_fac,
                                const RatesParams& params)
{
    Dual log_scor = actual_log_screen_d(state, scn_fac, params);
    Dual result = dual_const(1.0);
    if (!(log_scor.val <= 0.0)) {
        constexpr Real log_scor_max = 300.0;
        log_scor = dual_min(log_scor, dual_const(log_scor_max));

        result = dual_exp(log_scor);
    }
    return result;
}

} // namespace mpl

#endif
