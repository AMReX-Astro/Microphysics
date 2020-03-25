#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_REAL.H>
#include <AMReX_Algorithm.H>
#include <network_properties.H>
#include <screen.H>
#include <screen_data.H>

using namespace amrex;
using namespace scrn;

void screening_init() {

}

void screening_finalize() {

}


void add_screening_factor(const int i,
                          const Real z1, const Real a1, const Real z2, const Real a2) {

  BL_ASSERT(i < NSCREEN);

  scn_facs[i].z1 = z1;
  scn_facs[i].z2 = z2;
  scn_facs[i].a1 = a1;
  scn_facs[i].a2 = a2;

  scn_facs[i].zs13 = std::pow(scn_facs[i].z1 + scn_facs[i].z2, 1.0_rt/3.0_rt);
  scn_facs[i].zs13inv = 1.0_rt/scn_facs[i].zs13;
  scn_facs[i].zhat = std::pow(scn_facs[i].z1 + scn_facs[i].z2, 5.0_rt/3.0_rt) -
                     std::pow(scn_facs[i].z1, 5.0_rt/3.0_rt) - std::pow(scn_facs[i].z2, 5.0_rt/3.0_rt);
  scn_facs[i].zhat2 = std::pow(scn_facs[i].z1 + scn_facs[i].z2, 5.0_rt/12.0_rt) -
                     std::pow(scn_facs[i].z1, 5.0_rt/12.0_rt) - std::pow(scn_facs[i].z2, 5.0_rt/12.0_rt);
  scn_facs[i].lzav = (5.0_rt/3.0_rt) * std::log(scn_facs[i].z1 * scn_facs[i].z2 / (scn_facs[i].z1 + scn_facs[i].z2));
  scn_facs[i].aznut = std::pow(scn_facs[i].z1 * scn_facs[i].z1 *
                              scn_facs[i].z2 * scn_facs[i].z2 *
                              scn_facs[i].a1 * scn_facs[i].a2 /
                              (scn_facs[i].a1 + scn_facs[i].a2), 1.0_rt/3.0_rt);
}

AMREX_GPU_HOST_DEVICE
void screen5(const plasma_state_t state,
             const int jscreen,
             Real& scor, Real& scordt, Real& scordd) {

  // this subroutine calculates screening factors and their derivatives
  // for nuclear reaction rates in the weak, intermediate and strong regimes.
  // based on graboske, dewit, grossman and cooper apj 181 457 1973 for
  // weak screening. based on alastuey and jancovici apj 226 1034 1978,
  // with plasma parameters from itoh et al apj 234 1079 1979, for strong
  // screening.

  // input:
  // state   = plasma state (T, rho, abar, zbar, etc.)
  // jscreen = counter of which reaction is being calculated

  // output:
  // scor    = screening correction
  // scordt  = derivative of screening correction with temperature
  // scordd  = derivative of screening correction with density

  // Get the ion data based on the input index
  Real z1 = scn_facs[jscreen].z1;
  Real a1 = scn_facs[jscreen].a1;
  Real z2 = scn_facs[jscreen].z2;
  Real a2 = scn_facs[jscreen].a2;

  // calculate individual screening factors
  Real bb = z1 * z2;
  Real gamp = state.aa;
  Real gampdt = state.daadt;
  // Real gampdd = state.daadd;

  Real qq = fact * bb * scn_facs[jscreen].zs13inv;
  Real gamef = qq * gamp;
  Real gamefdt = qq * gampdt;
  // Real gamefdd  = qq * gampdd;

  Real tau12 = state.taufac * scn_facs[jscreen].aznut;
  Real tau12dt = state.taufacdt * scn_facs[jscreen].aznut;

  qq = 1.0_rt/tau12;
  Real alph12 = gamef * qq;
  Real alph12dt = (gamefdt - alph12*tau12dt) * qq;
  // Real alph12dd = gamefdd * qq;


  // limit alph12 to 1.6 to prevent unphysical behavior.
  // this should really be replaced by a pycnonuclear reaction rate formula
  if (alph12 > 1.6_rt) {
    alph12   = 1.6e0_rt;
    alph12dt = 0.0_rt;
    // alph12dd = 0.0_rt;

    gamef    = 1.6e0_rt * tau12;
    gamefdt  = 1.6e0_rt * tau12dt;
    // gamefdd  = 0.0_rt;

    qq = scn_facs[jscreen].zs13/(fact * bb);
    gamp = gamef * qq;
    gampdt = gamefdt * qq;
    // gampdd   = 0.0_rt;
  }

  // weak screening regime
  Real h12w = bb * state.qlam0z;
  Real dh12wdt = bb * state.qlam0zdt;
  //Real dh12wdd = bb * qlam0zdd;

  Real h12 = h12w;
  Real dh12dt = dh12wdt;
  //Real dh12dd  = dh12wdd;

  // intermediate and strong sceening regime
  if (gamef > gamefx) {

    Real gamp14 = std::pow(gamp, 0.25_rt);
    Real rr = 1.0_rt/gamp;
    qq = 0.25_rt * gamp14 * rr;
    Real gamp14dt = qq * gampdt;
    //Real gamp14dd = qq * gampdd;

    Real cc = 0.896434e0_rt * gamp * scn_facs[jscreen].zhat
      - 3.44740e0_rt * gamp14 * scn_facs[jscreen].zhat2
      - 0.5551e0_rt * (std::log(gamp) + scn_facs[jscreen].lzav)
      - 2.996e0_rt;

    Real dccdt = 0.896434e0_rt * gampdt * scn_facs[jscreen].zhat
      - 3.44740e0_rt * gamp14dt * scn_facs[jscreen].zhat2
      - 0.5551e0_rt *rr * gampdt;

    //dccdd    =   0.896434e0_rt * gampdd * zhat(jscreen) &
    //     - 3.44740e0_rt  * gamp14dd * zhat2(jscreen) &
    //     - 0.5551e0_rt*rr*gampdd

    Real a3 = alph12 * alph12 * alph12;
    Real da3 = 3.0e0_rt * alph12 * alph12;

    qq = 0.014e0_rt + 0.0128e0_rt*alph12;
    Real dqqdt  = 0.0128e0_rt*alph12dt;
    //dqqdd  = 0.0128e0_rt*alph12dd

    rr = (5.0_rt/32.0_rt) - alph12*qq;
    Real drrdt  = -(alph12dt*qq + alph12*dqqdt);
    // drrdd  = -(alph12dd*qq + alph12*dqqdd)

    Real ss = tau12*rr;
    Real dssdt  = tau12dt*rr + tau12*drrdt;
    // dssdd  = tau12*drrdd

    Real tt = -0.0098e0_rt + 0.0048e0_rt*alph12;
    Real dttdt  = 0.0048e0_rt*alph12dt;
    // dttdd  = 0.0048e0_rt*alph12dd

    Real uu = 0.0055e0_rt + alph12*tt;
    Real duudt  = alph12dt*tt + alph12*dttdt;
    // duudd  = alph12dd*tt + alph12*dttdd

    Real vv = gamef * alph12 * uu;
    Real dvvdt = gamefdt*alph12*uu + gamef*alph12dt*uu + gamef*alph12*duudt;
    // dvvdd= gamefdd*alph12*uu + gamef*alph12dd*uu + gamef*alph12*duudd

    h12 = cc - a3 * (ss + vv);
    rr = da3 * (ss + vv);
    dh12dt  = dccdt - rr*alph12dt - a3*(dssdt + dvvdt);
    // dh12dd  = dccdd - rr*alph12dd - a3*(dssdd + dvvdd)

    rr = 1.0_rt - 0.0562e0_rt*a3;
    ss = -0.0562e0_rt*da3;
    drrdt = ss*alph12dt;
    // drrdd  = ss*alph12dd

    Real xlgfac;
    Real dxlgfacdt;

    if (rr >= 0.77e0_rt) {
      xlgfac = rr;
      dxlgfacdt = drrdt;
      //dxlgfacdd = drrdd;
    } else {
      xlgfac = 0.77e0_rt;
      dxlgfacdt = 0.0_rt;
      //dxlgfacdd = 0.0_rt
    }

    h12 = std::log(xlgfac) + h12;
    rr = 1.0_rt/xlgfac;
    dh12dt = rr*dxlgfacdt + dh12dt;
    // dh12dd = rr*dxlgfacdd + dh12dd

    if (gamef <= gamefs) {
      Real dgamma = 1.0e0_rt/(gamefs - gamefx);

      rr =  dgamma*(gamefs - gamef);
      drrdt  = -dgamma*gamefdt;
      //drrdd  = -dgamma*gamefdd

      ss = dgamma*(gamef - gamefx);
      dssdt = dgamma*gamefdt;
      //dssdd  = dgamma*gamefdd

      vv = h12;

      h12 = h12w*rr + vv*ss;
      dh12dt = dh12wdt*rr + h12w*drrdt + dh12dt*ss + vv*dssdt;
      //dh12dd = dh12wdd*rr + h12w*drrdd + dh12dd*ss + vv*dssdd
    }

    // end of intermediate and strong screening
  }

  // machine limit the output
  // further limit to avoid the pycnonuclear regime
  h12 = amrex::max(amrex::min(h12, h12_max), 0.0_rt);
  scor = std::exp(h12);

  if (h12 == h12_max) {
    scordt = 0.0;
    //scordd = 0.0_rt
  } else {
    scordt = scor * dh12dt;
    //scordd = scor * dh12dd
  }
}


