#ifndef _screen_H_
#define _screen_H_

#include <AMReX.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_REAL.H>
#include <network_properties.H>
#include <microphysics_math.H>
#include <fundamental_constants.H>
#include <cmath>
#include <screen_data.H>
#include <extern_parameters.H>

using namespace amrex;

#if SCREEN_METHOD == 0
const std::string screen_name = "screen5";
#elif SCREEN_METHOD == 1
const std::string screen_name = "chugunov2007";
#elif SCREEN_METHOD == 2
const std::string screen_name = "chugunov2009";
#elif SCREEN_METHOD == 3
const std::string screen_name = "null";
#endif

struct plasma_state_t {

  Real qlam0z;
  Real qlam0zdt;
  //Real qlam0zdd;

  Real taufac;
  Real taufacdt;

  Real aa;
  Real daadt;
  //Real daadd;

  Real temp;
  Real zbar;
  Real z2bar;
  Real n_e;
  Real gamma_e_fac;
};

inline
std::ostream& operator<< (std::ostream& o, plasma_state_t const& pstate)
{
  o << "qlam0z   = " << pstate.qlam0z << std::endl;
  o << "qlam0zdt = " << pstate.qlam0zdt << std::endl;
  o << "taufac   = " << pstate.taufac << std::endl;
  o << "taufacdt = " << pstate.taufacdt << std::endl;
  o << "aa       = " << pstate.aa << std::endl;
  o << "daadt    = " << pstate.daadt << std::endl;
  o << "temp     = " << pstate.temp << std::endl;
  o << "zbar     = " << pstate.zbar << std::endl;
  o << "z2bar    = " << pstate.z2bar << std::endl;
  o << "n_e      = " << pstate.n_e << std::endl;
  o << "gamma_e_fac = " << pstate.gamma_e_fac << std::endl;

  return o;
}

AMREX_FORCE_INLINE
void
screening_init() {

}

AMREX_FORCE_INLINE
void
screening_finalize() {

}


AMREX_FORCE_INLINE
void add_screening_factor(const int i,
                          const Real z1, const Real a1, const Real z2, const Real a2)
{
  using namespace scrn;

  BL_ASSERT(i < NSCREEN);

  scn_facs[i].z1 = z1;
  scn_facs[i].z2 = z2;
  scn_facs[i].a1 = a1;
  scn_facs[i].a2 = a2;

  scn_facs[i].zs13 = std::cbrt(z1 + z2);
  scn_facs[i].zs13inv = 1.0_rt/scn_facs[i].zs13;
  scn_facs[i].zhat = std::pow(z1 + z2, 5.0_rt/3.0_rt) -
                     std::pow(z1, 5.0_rt/3.0_rt) -
                     std::pow(z2, 5.0_rt/3.0_rt);
  scn_facs[i].zhat2 = std::pow(z1 + z2, 5.0_rt/12.0_rt) -
                      std::pow(z1, 5.0_rt/12.0_rt) -
                      std::pow(z2, 5.0_rt/12.0_rt);
  scn_facs[i].lzav = (5.0_rt/3.0_rt) * std::log(z1 * z2 / (z1 + z2));
  scn_facs[i].aznut = std::cbrt(z1 * z1 * z2 * z2 * a1 * a2 / (a1 + a2));
  scn_facs[i].ztilde = 0.5_rt * (std::cbrt(z1) + std::cbrt(z2));

}


template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
fill_plasma_state(plasma_state_t& state, const Real temp, const Real dens, Array1D<Real, 1, NumSpec> const& y) {

  Real sum = 0.0_rt;
  for (int n = 1; n <= NumSpec; n++) {
    sum += y(n);
  }
  Real abar = 1.0_rt / sum;
  Real ytot = sum;

  sum = 0.0_rt;
  Real sum2 = 0.0_rt;
  for (int n = 1; n <= NumSpec; n++) {
    sum += zion[n-1]*y(n);
    sum2 += zion[n-1]*zion[n-1]*y(n);
  }

  const Real co2 = (1.0_rt/3.0_rt) * 4.248719e3_rt;

  Real zbar = sum * abar;
  Real z2bar = sum2 * abar;

  // ntot
  Real rr = dens * ytot;
  Real tempi = 1.0_rt / temp;
  Real dtempi;
  if constexpr (do_T_derivatives) {
    dtempi = -tempi * tempi;
  }

  Real pp = std::sqrt(rr*tempi*(z2bar + zbar));
  Real dppdt;
  if constexpr (do_T_derivatives) {
    Real qq = 0.5_rt/pp *(z2bar + zbar);
    dppdt = qq*rr*dtempi;
  }

  state.qlam0z = 1.88e8_rt * tempi * pp;
  if constexpr (do_T_derivatives) {
    state.qlam0zdt = 1.88e8_rt * (dtempi*pp + tempi*dppdt);
  }

  state.taufac = co2 * std::cbrt(tempi);
  if constexpr (do_T_derivatives) {
    state.taufacdt = -(1.0_rt/3.0_rt) * state.taufac * tempi;
  }

  Real xni = std::cbrt(rr * zbar);

  state.aa = 2.27493e5_rt * tempi * xni;
  if constexpr (do_T_derivatives) {
    state.daadt = 2.27493e5_rt * dtempi * xni;
  }

  state.temp = temp;
  state.zbar = zbar;
  state.z2bar = z2bar;

  // Electron number density
  // zbar * ntot works out to sum(z[i] * n[i]), after cancelling terms
  state.n_e = zbar * rr * C::n_A;

  // precomputed part of Gamma_e, from Chugunov 2009 eq. 6
  constexpr Real gamma_e_constants =
      C::q_e*C::q_e/C::k_B * gcem::pow(4.0_rt/3.0_rt*M_PI, 1.0_rt/3.0_rt);
  state.gamma_e_fac = gamma_e_constants * std::cbrt(state.n_e);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
fill_plasma_state(plasma_state_t& state, const Real temp, const Real dens, Array1D<Real, 1, NumSpec> const& y) {
    if (jacobian == 1) {
        constexpr int do_T_derivatives = 1;
        fill_plasma_state<do_T_derivatives>(state, temp, dens, y);
    } else {
        constexpr int do_T_derivatives = 0;
        fill_plasma_state<do_T_derivatives>(state, temp, dens, y);
    }
}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_screen5 (const plasma_state_t& state,
                     const scrn::screen_factors_t& scn_fac,
                     Real& scor, Real& scordt)
{
    // this subroutine calculates screening factors and their derivatives
    // for nuclear reaction rates in the weak, intermediate and strong regimes.
    // based on graboske, dewit, grossman and cooper apj 181 457 1973 for
    // weak screening. based on alastuey and jancovici apj 226 1034 1978,
    // with plasma parameters from itoh et al apj 234 1079 1979, for strong
    // screening.

    // input:
    // state   = plasma state (T, rho, abar, zbar, etc.)
    // scn_fac = screening factors for A and Z

    // output:
    // scor    = screening correction
    // scordt  = derivative of screening correction with temperature

    const Real fact    = 1.25992104989487e0_rt;
    const Real gamefx  = 0.3e0_rt;
    const Real gamefs  = 0.8e0_rt;
    const Real h12_max = 300.e0_rt;

    // Get the ion data based on the input index
    Real z1 = scn_fac.z1;
    Real z2 = scn_fac.z2;

    // calculate individual screening factors
    Real bb = z1 * z2;
    Real gamp = state.aa;
    Real gampdt;
    if constexpr (do_T_derivatives) {
        gampdt = state.daadt;
    }

    Real qq = fact * bb * scn_fac.zs13inv;
    Real gamef = qq * gamp;
    Real gamefdt;
    if constexpr (do_T_derivatives) {
        gamefdt = qq * gampdt;
    }

    Real tau12 = state.taufac * scn_fac.aznut;
    Real tau12dt;
    if constexpr (do_T_derivatives) {
        tau12dt = state.taufacdt * scn_fac.aznut;
    }

    qq = 1.0_rt/tau12;
    Real alph12 = gamef * qq;
    Real alph12dt;
    if constexpr (do_T_derivatives) {
        alph12dt = (gamefdt - alph12*tau12dt) * qq;
    }


    // limit alph12 to 1.6 to prevent unphysical behavior.
    // this should really be replaced by a pycnonuclear reaction rate formula
    if (alph12 > 1.6_rt) {
        alph12   = 1.6e0_rt;
        if constexpr (do_T_derivatives) {
            alph12dt = 0.0_rt;
        }

        gamef    = 1.6e0_rt * tau12;
        if constexpr (do_T_derivatives) {
            gamefdt  = 1.6e0_rt * tau12dt;
        }

        qq = scn_fac.zs13/(fact * bb);
        gamp = gamef * qq;
        if constexpr (do_T_derivatives) {
            gampdt = gamefdt * qq;
        }
    }

    // weak screening regime
    Real h12w = bb * state.qlam0z;
    Real dh12wdt;
    if constexpr (do_T_derivatives) {
        dh12wdt = bb * state.qlam0zdt;
    }

    Real h12 = h12w;
    Real dh12dt;
    if constexpr (do_T_derivatives) {
        dh12dt = dh12wdt;
    }

    // intermediate and strong sceening regime
    if (gamef > gamefx) {

        Real gamp14 = std::pow(gamp, 0.25_rt);
        Real rr = 1.0_rt/gamp;

        Real cc = 0.896434e0_rt * gamp * scn_fac.zhat
            - 3.44740e0_rt * gamp14 * scn_fac.zhat2
            - 0.5551e0_rt * (std::log(gamp) + scn_fac.lzav)
            - 2.996e0_rt;

        Real dccdt;
        if constexpr (do_T_derivatives) {
            qq = 0.25_rt * gamp14 * rr;
            Real gamp14dt = qq * gampdt;
            dccdt = 0.896434e0_rt * gampdt * scn_fac.zhat
                - 3.44740e0_rt * gamp14dt * scn_fac.zhat2
                - 0.5551e0_rt *rr * gampdt;
        }

        Real a3 = alph12 * alph12 * alph12;
        Real da3 = 3.0e0_rt * alph12 * alph12;

        qq = 0.014e0_rt + 0.0128e0_rt*alph12;

        rr = (5.0_rt/32.0_rt) - alph12*qq;
        Real drrdt;
        if constexpr (do_T_derivatives) {
            Real dqqdt  = 0.0128e0_rt*alph12dt;
            drrdt  = -(alph12dt*qq + alph12*dqqdt);
        }

        Real ss = tau12*rr;

        Real tt = -0.0098e0_rt + 0.0048e0_rt*alph12;

        Real uu = 0.0055e0_rt + alph12*tt;

        Real vv = gamef * alph12 * uu;

        h12 = cc - a3 * (ss + vv);
        if constexpr (do_T_derivatives) {
            Real dssdt  = tau12dt*rr + tau12*drrdt;
            Real dttdt  = 0.0048e0_rt*alph12dt;
            Real duudt  = alph12dt*tt + alph12*dttdt;
            Real dvvdt = gamefdt*alph12*uu + gamef*alph12dt*uu + gamef*alph12*duudt;
            rr = da3 * (ss + vv);
            dh12dt  = dccdt - rr*alph12dt - a3*(dssdt + dvvdt);
        }

        rr = 1.0_rt - 0.0562e0_rt*a3;
        if constexpr (do_T_derivatives) {
            ss = -0.0562e0_rt*da3;
            drrdt = ss*alph12dt;
        }

        Real xlgfac;
        Real dxlgfacdt;

        if (rr >= 0.77e0_rt) {
            xlgfac = rr;
            if constexpr (do_T_derivatives) {
                dxlgfacdt = drrdt;
            }
        } else {
            xlgfac = 0.77e0_rt;
            if constexpr (do_T_derivatives) {
                dxlgfacdt = 0.0_rt;
            }
        }

        h12 = std::log(xlgfac) + h12;
        rr = 1.0_rt/xlgfac;
        if constexpr (do_T_derivatives) {
            dh12dt = rr*dxlgfacdt + dh12dt;
        }

        if (gamef <= gamefs) {
            Real dgamma = 1.0e0_rt/(gamefs - gamefx);

            rr =  dgamma*(gamefs - gamef);

            ss = dgamma*(gamef - gamefx);

            vv = h12;

            h12 = h12w*rr + vv*ss;
            if constexpr (do_T_derivatives) {
                drrdt  = -dgamma*gamefdt;
                Real dssdt = dgamma*gamefdt;
                dh12dt = dh12wdt*rr + h12w*drrdt + dh12dt*ss + vv*dssdt;
            }
        }

        // end of intermediate and strong screening
    }

    // machine limit the output
    // further limit to avoid the pycnonuclear regime
    h12 = amrex::max(amrex::min(h12, h12_max), 0.0_rt);
    scor = std::exp(h12);

    if constexpr (do_T_derivatives) {
        if (h12 == h12_max) {
            scordt = 0.0;
        } else {
            scordt = scor * dh12dt;
        }
    }
}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void chugunov2007 (const plasma_state_t& state,
                   const scrn::screen_factors_t& scn_fac,
                   Real& scor, Real& scordt)
{
    // Calculates screening factors based on Chugunov et al. 2007, following the
    // the approach in Yakovlev 2006 to extend to a multi-component plasma.
    //
    // References:
    //   Chugunov, DeWitt, and Yakovlev 2007, PhRvD, 76, 025028
    //   Yakovlev, Gasques, Afanasjev, Beard, and Wiescher 2006, PhRvC, 74, 035803
    //   Chugunov and DeWitt 2009, PhRvC, 80, 014611

    // input:
    // state   = plasma state (T, rho, abar, zbar, etc.)
    // scn_fac = screening factors for A and Z

    // output:
    // scor    = screening correction
    // scordt  = derivative of screening correction with temperature

    Real tmp;
    // Plasma temperature T_p
    // This formula comes from working backwards from zeta_ij (Chugunov 2009 eq. 12)
    // through Chugunov 2007 eq. 3 to Chugunov 2007 eq. 2.
    // Ultimately, the changes from the expression in Chugunov 2007 are:
    //   Z^2 -> Z1 * Z2
    //   n_i -> n_e / ztilde^3, where ztilde = (Z1^(1/3) + Z2^(1/3)) / 2
    //   m_i -> 2 mu12 (reduced mass)
    // This prescription reduces to the expressions from Chugunov 2007 in the case
    // of an OCP, and to Chugunov 2009 in the case of a binary ionic mixture.
    // This also matches Yakovlev et al. 2006, eq. 10.
    //
    // For reference, MESA r21.12.1 does:
    //   Z^2 -> Z1 * Z2
    //   n_i -> n_e / zbar (=ntot)
    //   m_i -> m_u * abar
    // Sam Jones' Fortran implementation (linked in issue #146) does:
    //   Z^2 -> zbar^2
    //   n_i -> ntot
    //   m_i -> m_u * abar
    Real mu12 = scn_fac.a1 * scn_fac.a2 / (scn_fac.a1 + scn_fac.a2);
    Real z_factor = scn_fac.z1 * scn_fac.z2;
    Real n_i = state.n_e / gcem::pow(scn_fac.ztilde, 3);
    Real m_i = 2.0_rt * mu12 / C::n_A;

    constexpr Real T_p_factor = C::hbar/C::k_B*C::q_e*gcem::sqrt(4.0_rt*GCEM_PI);
    Real T_p = T_p_factor * std::sqrt(z_factor * n_i / m_i);

    // Normalized temperature
    Real inv_T_p = 1.0_rt / T_p;
    Real T_norm = state.temp * inv_T_p;
    Real dT_norm_dT = inv_T_p;

    // The fit has only been verified down to T ~ 0.1 T_p, below which the rate
    // should be nearly temperature-independent (in the pycnonuclear regime),
    // and we clip the temperature to 0.1 T_p at small T.
    // start the transition here
    constexpr Real T_norm_fade = 0.2_rt;
    constexpr Real T_norm_min = 0.1_rt;

    if (T_norm < T_norm_min) {
        // clip temperature to the minimum value
        T_norm = T_norm_min;
        dT_norm_dT = 0.0_rt;
    } else if (T_norm <= T_norm_fade) {
        // blend using a cosine, after MESA
        constexpr Real delta_T = T_norm_fade - T_norm_min;
        tmp = M_PI * (T_norm - T_norm_min) / delta_T;
        Real f = 0.5_rt * (1.0_rt - std::cos(tmp));
        Real df_dT;
        if constexpr (do_T_derivatives) {
            df_dT = 0.5_rt * M_PI / delta_T * std::sin(tmp) * dT_norm_dT;
        }

        if constexpr (do_T_derivatives) {
            dT_norm_dT = -df_dT * T_norm_min + df_dT * T_norm + f * dT_norm_dT;
        }
        T_norm = (1.0_rt - f) * T_norm_min + f * T_norm;
    }

    // Coulomb coupling parameter from Yakovlev 2006 eq. 10
    Real Gamma = state.gamma_e_fac*scn_fac.z1*scn_fac.z2 / (scn_fac.ztilde*T_norm*T_p);
    Real dGamma_dT;
    if constexpr (do_T_derivatives) {
        dGamma_dT = -Gamma / T_norm * dT_norm_dT;
    }

    // The fit for Gamma is only applicable up to ~600, so smoothly cap its value
    constexpr Real Gamma_fade = 590;
    constexpr Real Gamma_max = 600;

    if (Gamma > Gamma_max) {
        // clip Gamma to the max value
        Gamma = Gamma_max;
        dGamma_dT = 0.0_rt;
    } else if (Gamma >= Gamma_fade) {
        // blend using a cosine, after MESA
        constexpr Real delta_gamma = Gamma_max - Gamma_fade;
        tmp = M_PI * (Gamma - Gamma_fade) / delta_gamma;
        Real f = 0.5_rt * (1.0_rt - std::cos(tmp));
        Real df_dT;
        if constexpr (do_T_derivatives) {
            df_dT = 0.5_rt * M_PI / delta_gamma * std::sin(tmp) * dGamma_dT;
            dGamma_dT = dGamma_dT - (df_dT * Gamma + f * dGamma_dT) + df_dT * Gamma_max;
        }
        Gamma = (1.0_rt - f) * Gamma + f * Gamma_max;
    }

    // Chugunov 2007 eq. 3
    Real zeta = std::cbrt(4.0_rt / (3.0_rt * GCEM_PI*GCEM_PI * T_norm*T_norm));
    Real dzeta_dT;
    if constexpr (do_T_derivatives) {
        dzeta_dT = -2.0_rt / (3.0_rt * T_norm) * zeta * dT_norm_dT;
    }

    // Gamma tilde from Chugunov 2007 eq. 21
    constexpr Real fit_alpha = 0.022_rt;
    Real fit_beta = 0.41_rt - 0.6_rt / Gamma;
    Real fit_gamma = 0.06_rt + 2.2_rt / Gamma;
    // Polynomial term in Gamma tilde
    Real poly = 1.0_rt + zeta*(fit_alpha + zeta*(fit_beta + fit_gamma*zeta));
    Real dpoly_dT;
    if constexpr (do_T_derivatives) {
        tmp = dGamma_dT / (Gamma * Gamma);
        Real dfit_beta_dT = 0.6_rt * tmp;
        Real dfit_gamma_dT = -2.2_rt * tmp;
        dpoly_dT = (fit_alpha + 2.0_rt*zeta*fit_beta + 3.0_rt*fit_gamma*zeta*zeta) * dzeta_dT
                   + zeta*zeta*(dfit_beta_dT + dfit_gamma_dT*zeta);
    }

    Real gamtilde = Gamma / std::cbrt(poly);
    // this is gamtilde * dlog(gamtilde)/dT
    Real dgamtilde_dT;
    if constexpr (do_T_derivatives) {
        dgamtilde_dT = gamtilde * (dGamma_dT / Gamma - dpoly_dT / poly / 3.0_rt);
    }

    // fit parameters just after Chugunov 2007 eq. 19
    constexpr Real A1 = 2.7822_rt;
    constexpr Real A2 = 98.34_rt;
    constexpr Real A3 = gcem::sqrt(3.0_rt) - A1 / gcem::sqrt(A2);
    const Real B1 = -1.7476_rt;
    const Real B2 = 66.07_rt;
    const Real B3 = 1.12_rt;
    const Real B4 = 65_rt;
    Real gamtilde2 = gamtilde * gamtilde;

    Real term1, term2, term3, term4;
    Real dterm1_dT, dterm2_dT, dterm3_dT, dterm4_dT;
    // Chugunov 2007 eq. 19
    term1 = 1.0_rt / std::sqrt(A2 + gamtilde);
    if constexpr (do_T_derivatives) {
        dterm1_dT = -0.5_rt * term1 / (A2 + gamtilde) * dgamtilde_dT;
    }

    term2 = 1.0_rt / (1.0_rt + gamtilde);
    if constexpr (do_T_derivatives) {
        dterm2_dT = -term2 * term2 * dgamtilde_dT;
    }

    tmp = B2 + gamtilde;
    term3 = gamtilde2 / tmp;
    if constexpr (do_T_derivatives) {
        dterm3_dT = gamtilde * (B2 + tmp) / (tmp * tmp) * dgamtilde_dT;
    }

    tmp = B4 + gamtilde2;
    term4 = gamtilde2 / tmp;
    if constexpr (do_T_derivatives) {
        dterm4_dT = B4 / (tmp * tmp) * 2.0_rt * gamtilde * dgamtilde_dT;
    }

    Real inner = A1 * term1 + A3 * term2;

    Real gamtilde32 = std::pow(gamtilde, 1.5_rt);
    Real h = gamtilde32 * inner + B1 * term3 + B3 * term4;
    Real dh_dT;
    if constexpr (do_T_derivatives) {
        Real dinner_dT = A1 * dterm1_dT + A3 * dterm2_dT;
        Real dgamtilde32_dT = 1.5_rt * std::sqrt(gamtilde) * dgamtilde_dT;
        dh_dT = dgamtilde32_dT * inner + gamtilde32 * dinner_dT
                + B1 * dterm3_dT + B3 * dterm4_dT;
    }

    // machine limit the output
    constexpr Real h_max = 300.e0_rt;
    h = amrex::min(h, h_max);
    scor = std::exp(h);

    if constexpr (do_T_derivatives) {
        if (h == h_max) {
            scordt = 0.0;
        } else {
            scordt = scor * dh_dT;
        }
    }
}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void f0 (const Real gamma, const Real dlog_dT, Real& f, Real& df_dT)
{
    // Calculates the free energy per ion in a OCP, from Chugunov and DeWitt 2009
    // equation 24.

    // input:
    // gamma
    // dlog_dT = derivative of log(gamma) with temperature

    // output:
    // f       = free energy
    // df_dT   = derivative of free energy with temperature

    // fit parameters
    constexpr Real A1 = -0.907_rt;
    constexpr Real A2 = 0.62954_rt;
    constexpr Real A3 = -gcem::sqrt(3.0_rt) / 2.0_rt - A1 / gcem::sqrt(A2);
    constexpr Real B1 = 0.00456_rt;
    constexpr Real B2 = 211.6_rt;
    constexpr Real B3 = -1e-4_rt;
    constexpr Real B4 = 0.00462_rt;
    Real gamma_12 = std::sqrt(gamma);

    Real term1, term2, term3, term4, term5;
    Real dterm1_dgamma, dterm2_dgamma, dterm3_dgamma, dterm4_dgamma, dterm5_dgamma;

    term1 = gamma_12 * std::sqrt(A2 + gamma);
    if constexpr (do_T_derivatives) {
        dterm1_dgamma = (0.5_rt * A2 + gamma) / term1;
    }

    term2 = std::log(std::sqrt(gamma / A2) + std::sqrt(1.0_rt + gamma / A2));
    if constexpr (do_T_derivatives) {
        dterm2_dgamma = 1.0_rt / (2.0_rt * term1);
    }

    term3 = gamma_12 - std::atan(gamma_12);
    if constexpr (do_T_derivatives) {
        dterm3_dgamma = 0.5_rt * gamma_12 / (1.0_rt + gamma);
    }

    term4 = std::log1p(gamma / B2);
    if constexpr (do_T_derivatives) {
        dterm4_dgamma = 1.0_rt / (B2 + gamma);
    }

    term5 = std::log1p(gamma * gamma / B4);
    if constexpr (do_T_derivatives) {
        dterm5_dgamma = 2.0_rt * gamma / (B4 + gamma * gamma);
    }

    f = A1*(term1 - A2*term2) + 2.0_rt*A3*term3 + B1*(gamma - B2*term4) + 0.5_rt*B3*term5;
    if constexpr (do_T_derivatives) {
        df_dT = (A1*(dterm1_dgamma - A2*dterm2_dgamma) +
                 2.0_rt*A3*dterm3_dgamma +
                 B1*(1.0_rt - B2*dterm4_dgamma) +
                 0.5_rt*B3*dterm5_dgamma) * dlog_dT * gamma;
    }
}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void chugunov2009 (const plasma_state_t& state,
                   const scrn::screen_factors_t& scn_fac,
                   Real& scor, Real& scordt)
{
    // Calculates screening factors based on Chugunov and DeWitt 2009, PhRvC, 80, 014611

    // input:
    // state   = plasma state (T, rho, abar, zbar, etc.)
    // scn_fac = screening factors for A and Z

    // output:
    // scor    = screening correction
    // scordt  = derivative of screening correction with temperature

    Real z1z2 = scn_fac.z1 * scn_fac.z2;
    Real zcomp = scn_fac.z1 + scn_fac.z2;

    // Gamma_e from eq. 6
    Real Gamma_e = state.gamma_e_fac / state.temp;
    // work in terms of log derivatives, since it's much simpler
    Real dlog_Gamma_dT = -1.0_rt / state.temp;

    // Coulomb coupling parameters for ions and compound nucleus, eqs. 7 & 9
    Real Gamma_1 = Gamma_e * std::pow(scn_fac.z1, 5.0_rt/3.0_rt);
    Real Gamma_2 = Gamma_e * std::pow(scn_fac.z2, 5.0_rt/3.0_rt);
    Real Gamma_comp = Gamma_e * std::pow(zcomp, 5.0_rt/3.0_rt);

    Real Gamma_12 = Gamma_e * z1z2 / scn_fac.ztilde;

    // Coulomb barrier penetrability, eq. 10
    constexpr Real tau_factor = gcem::pow(
            27.0_rt/2.0_rt * gcem::pow(M_PI*C::q_e*C::q_e/C::hbar, 2) / (C::n_A*C::k_B),
            1.0_rt/3.0_rt);
    Real tau_12 = tau_factor * scn_fac.aznut / std::cbrt(state.temp);
    Real dlog_tau_12_dT;
    if constexpr (do_T_derivatives) {
        dlog_tau_12_dT = -1.0_rt / state.temp / 3.0_rt;
    }

    // eq. 12
    Real zeta = 3.0_rt * Gamma_12 / tau_12;
    Real dzeta_dT;
    if constexpr (do_T_derivatives) {
        dzeta_dT = zeta * (dlog_Gamma_dT - dlog_tau_12_dT);
    }

    // additional fit parameters, eq. 25
    Real y_12 = 4.0_rt * z1z2 / (zcomp * zcomp);
    Real c1 = 0.013_rt * y_12 * y_12;
    Real c2 = 0.406_rt * std::pow(y_12, 0.14_rt);
    Real c3 = 0.062_rt * std::pow(y_12, 0.19_rt) + 1.8_rt / Gamma_12;

    Real poly = 1.0_rt + zeta*(c1 + zeta*(c2 + c3*zeta));
    Real t_12 = std::cbrt(poly);
    Real dlog_dT;
    if constexpr (do_T_derivatives) {
        Real dc3_dT = -1.8_rt / Gamma_12 * dlog_Gamma_dT;
        Real dpoly_dT = (c1 + zeta*(2.0_rt*c2 + 3.0_rt*c3*zeta))*dzeta_dT + dc3_dT*zeta*zeta*zeta;
        Real dlog_t_12_dT = dpoly_dT / (3.0_rt * poly);
        dlog_dT = dlog_Gamma_dT - dlog_t_12_dT;
    }

    // strong screening enhancement factor, eq. 23, replacing tau_ij with t_ij
    // Using Gamma/tau_ij gives extremely low values, while Gamma/t_ij gives
    // values similar to those from Chugunov 2007.
    Real term1, term2, term3;
    Real dterm1_dT, dterm2_dT, dterm3_dT;
    f0<do_T_derivatives>(Gamma_1 / t_12, dlog_dT, term1, dterm1_dT);
    f0<do_T_derivatives>(Gamma_2 / t_12, dlog_dT, term2, dterm2_dT);
    f0<do_T_derivatives>(Gamma_comp / t_12, dlog_dT, term3, dterm3_dT);
    Real h_fit = term1 + term2 - term3;
    Real dh_fit_dT;
    if constexpr (do_T_derivatives) {
        dh_fit_dT = dterm1_dT + dterm2_dT - dterm3_dT;
    }

    // weak screening correction term, eq. A3
    Real corr_C = 3.0_rt*z1z2 * std::sqrt(state.z2bar/state.zbar) /
                  (std::pow(zcomp, 2.5_rt) - std::pow(scn_fac.z1, 2.5_rt) - std::pow(scn_fac.z2, 2.5_rt));

    // corrected enhancement factor, eq. A4
    Real Gamma_12_2 = Gamma_12 * Gamma_12;
    Real numer = corr_C + Gamma_12_2;
    Real denom = 1.0_rt + Gamma_12_2;
    Real h12 = numer / denom * h_fit;
    Real dh12_dT;
    if constexpr (do_T_derivatives) {
        Real dGamma_12_2_dT = 2 * Gamma_12_2 * dlog_Gamma_dT;
        dh12_dT = h12 * (dGamma_12_2_dT/numer - dGamma_12_2_dT/denom + dh_fit_dT/h_fit);
    }

    // machine limit the output
    const Real h12_max = 300.e0_rt;
    h12 = amrex::min(h12, h12_max);
    scor = std::exp(h12);

    if constexpr (do_T_derivatives) {
        if (h12 == h12_max) {
            scordt = 0.0_rt;
        } else {
            scordt = scor * dh12_dT;
        }
    }
}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void actual_screen(const plasma_state_t& state,
                   const scrn::screen_factors_t& scn_fac,
                   Real& scor, Real& scordt)
{
#if SCREEN_METHOD == 0
    actual_screen5<do_T_derivatives>(state, scn_fac, scor, scordt);
#elif SCREEN_METHOD == 1
    chugunov2007<do_T_derivatives>(state, scn_fac, scor, scordt);
#elif SCREEN_METHOD == 2
    chugunov2009<do_T_derivatives>(state, scn_fac, scor, scordt);
#elif SCREEN_METHOD == 3
    // null screening
    scor = 1.0_rt;
    scordt = 0.0_rt;
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void actual_screen(const plasma_state_t& state,
                   const scrn::screen_factors_t& scn_fac,
                   Real& scor, Real& scordt)
{
    if (jacobian == 1) {
        constexpr int do_T_derivatives = 1;
        actual_screen<do_T_derivatives>(state, scn_fac, scor, scordt);
    } else {
        constexpr int do_T_derivatives = 0;
        actual_screen<do_T_derivatives>(state, scn_fac, scor, scordt);
        scordt = 0.0_rt;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void screen(const plasma_state_t& state,
            const int jscreen,
            MICROPHYSICS_UNUSED const Real z1_screen, MICROPHYSICS_UNUSED const Real a1_screen,
            MICROPHYSICS_UNUSED const Real z2_screen, MICROPHYSICS_UNUSED const Real a2_screen,
            Real& scor, Real& scordt, Real& /*scordd*/)
{
    using namespace scrn;

#ifndef AMREX_USE_GPU
    AMREX_ASSERT(scn_facs[jscreen].validate_nuclei(z1_screen, a1_screen, z2_screen, a2_screen));
#endif

    actual_screen(state, scn_facs[jscreen], scor, scordt);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void screen(const plasma_state_t& state,
            const scrn::screen_factors_t& scn_fac,
            Real& scor, Real& scordt, Real& /*scordd*/)
{
    using namespace scrn;

    actual_screen(state, scn_fac, scor, scordt);
}

#endif