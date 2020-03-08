#ifndef _actual_conductivity_H_
#define _actual_conductivity_H_

#include <cmath>
#include "eos_type.H"
#include "network.H"
#include "fundamental_constants.H"

AMREX_FORCE_INLINE
void
actual_conductivity_init()
{

}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
actual_conductivity(eos_t state) {

  // this routine is sig99, it approximates the thermal transport
  // coefficients.
  //
  // input:
  // temp   = temperature temp (in K) = state % T
  // den    = density den (in g/cm**3) = state % rho
  // nspec  = number of isotopes in the composition
  // xmass  = mass fractions of the composition = state % xn
  // zion   = number of protons in each isotope (charge of each isotope)
  // aion   = number of protons + neutrons in each isotope (atomic weight)
  // pep    = electron-positron pressure (in erg/cm**3) = state % pele
  // xne    = electron-positron number density (in 1/cm**3) = state % xne
  // eta    = electron degeneracy parameter (chemical potential / k T) = state % eta

  // output:
  // orad   = radiation contribution to the opacity (in cm**2/g)
  // ocond  = conductive contribution to the opacity (in cm**2/g)
  // opac   = the total opacity (in cm**2/g)
  // conductivity = thermal conductivity (in erg/cm/K/sec)
  // 

    // various physical and derived constants
    // con2 = con1*sqrt(4*pi*e*e/me)
    // meff = hbar/(me*c)*(3*pi**2)**(1/3)
    // weid = (pi*kerg)**2/(3*me)
    // iec  = 4*e**4*me/(3*pi*hbar**3)
    // xec  = hbar/kerg*e*sqrt(4*pi/me)

  const Real third  = 1.0_rt/3.0_rt;
  const Real twoth  = 2.0_rt * third;

  const Real zbound = 0.1e0_rt;
  const Real t7peek = 1.0e20_rt;
  const Real k2c    = 4.0_rt/3.0_rt*a_rad*c_light;
  const Real meff   = 1.194648642401440e-10_rt;
  const Real weid   = 6.884326138694269e-5_rt;
  const Real iec    = 1.754582332329132e16_rt;
  const Real xec    = 4.309054377592449e-7_rt;
  const Real rt3    = 1.7320508075688772e0_rt;
  const Real con2   = 1.07726359439811217e-7_rt;

  // switches for the iben & christy regimes
  const Real t6_switch1 = 0.5_rt;
  const Real t6_switch2 = 0.9_rt;

  // initialize
  Real opac      = 0.0e0_rt;
  Real orad      = 0.0e0_rt;
  Real ocond     = 0.0e0_rt;
  Real oiben1    = 0.0e0_rt;
  Real oiben2    = 0.0e0_rt;
  Real ochrs     = 0.0e0_rt;
  Real oh        = 0.0e0_rt;
  Real ov        = 0.0e0_rt;
  Real zbar      = 0.0e0_rt;
  Real ytot1     = 0.0e0_rt;

  // set the composition variables
  Real w[6];
  for (int i = 0; i < 6; i++) {
    w[i] = 0.0e0_rt;
  }

  // the idea here is that w[0] is H, w[1] is He, and w[2] is metals
  for (int i = 0; i < NumSpec; i++) {
    int iz = amrex::min(3, amrex::max(1, static_cast<int>(zion[i]))) - 1;
    Real ymass = state.xn[i]*aion_inv[i];
    w[iz] += state.xn[i];
    w[iz+3] += zion[i] * zion[i] * ymass;
    zbar += zion[i] * ymass;
    ytot1 += ymass;
  }
  Real abar = 1.0e0_rt/ytot1;
  Real zbar = zbar * abar;
  Real t6 = state.T * 1.0e-6_rt;

  Real xh = w[0];
  Real xhe = w[1];
  Real xz = w[2];

  // radiative section:
  // from iben apj 196 525 1975
  if (xh < 1.0e-5_rt) {
    Real xmu = amrex::max(1.0e-99_rt, w[3] + w[4] + w[5] - 1.0e0_rt);
    Real xkc = std::pow((2.019e-4_rt * state.rho / std::pow(t6, 1.7_rt)), 2.425_rt);
    Real xkap = 1.0_rt + xkc * (1.0_rt + xkc/24.55_rt);
    Real xkb = 3.86_rt + 0.252_rt*std::sqrt(xmu) + 0.018_rt*xmu;
    Real xka = 3.437_rt * (1.25_rt + 0.488_rt*std::sqrt(xmu) + 0.092_rt*xmu);
    Real dbar = std::exp(-xka + xkb*std::log(t6));
    oiben1 = xkap * std::pow(state.rho/dbar, 0.67_rt);
  }

  if ( !((xh >=  1.0e-5_rt) && (t6 < t6_switch1)) &&
       !((xh < 1.0e-5_rt) && (xz > zbound)) ) {
    Real d0log;
    if (t6 > t6_switch1) {
      d0log = -(3.868_rt + 0.806_rt*xh) + 1.8_rt*std::log(t6);
    } else {
      d0log = -(3.868_rt + 0.806_rt*xh) + (3.42_rt - 0.52_rt*xh)*std::log(t6);
    }
    Real xka1 = 2.809_rt * std::exp(-(1.74_rt  - 0.755_rt*xh)
                                    * std::pow(std::log10(t6) - 0.22_rt + 0.1375_rt*xh, 2));

   Real xkw = 4.05_rt * std::exp(-(0.306_rt  - 0.04125_rt*xh)
                                 * std::pow(std::log10(t6) - 0.18_rt + 0.1625_rt*xh, 2));
   Real xkaz = 50.0_rt*xz*xka1 * std::exp(-0.5206_rt*std::pow((std::log(state.rho)-d0log)/xkw, 2));
   Real dbar2log = -(4.283_rt + 0.7196_rt*xh) + 3.86_rt*std::log(t6);
   Real dbar1log = -5.296_rt + 4.833_rt*std::log(t6);
   if (dbar2log < dbar1log) {
     dbar1log = dbar2log;
   }
   oiben2 = std::pow(state.rho/std::exp(dbar1log), 0.67_rt) * std::exp(xkaz);
  }

  // from christy apj 144 108 1966
  if ((t6 < t6_switch2) && (xh >= 1.0e-5_rt)) {
    Real t4 = state % T * 1.0e-4_rt;
    Real t4r = std::sqrt(t4);
    Real t44 = t4*t4*t4*t4;
    Real t45 = t44 * t4;
    Real t46 = t45 * t4;
    Real ck1 = 2.0e6_rt/t44 + 2.1_rt*t46;
    Real ck3 = 4.0e-3_rt/t44 + 2.0e-4_rt/std::pow(state.rho, 0.25_rt);
    Real ck2 = 4.5_rt*t46 + 1.0_rt/(t4*ck3);
    Real ck4 = 1.4e3_rt*t4 + t46;
    Real ck5 = 1.0e6_rt + 0.1_rt*t46;
    Real ck6 = 20.0_rt*t4 + 5.0_rt*t44 + t45;
    Real xkcx = xh*(t4r/ck1 + 1.0_rt/ck2);
    Real xkcy = xhe*(1.0_rt/ck4 + 1.5_rt/ck5);
    Real xkcz = xz*(t4r/ck6);
    ochrs = state.pele * (xkcx + xkcy + xkcz);
  }

  // opacity in presence of hydrogen
  if (xh >= 1.0e-5_rt) {
    if (t6 < t6_switch1) {
      orad = ochrs;
    } else if (t6 <= t6_switch2) {
      orad = 2.0_rt*(ochrs*(1.5_rt - t6) + oiben2*(t6 - 1.0_rt));
    } else {
      orad = oiben2;
    }

    // opacity in absence of hydrogen
  } else {
    if (xz > zbound) {
      orad = oiben1;
    } else {
      orad = oiben1*(xz/zbound) + oiben2*((zbound-xz)/zbound);
    }
  }

  // add in the compton scattering opacity, weaver et al. apj 1978 225 1021
  Real th = std::min(511.0_rt, state.T * 8.617e-8_rt);
  Real fact = 1.0_rt + 2.75e-2_rt*th - 4.88e-5_rt*th*th;
  Real facetax = 1.0e100_rt;
  if (state.eta <= 500.0_rt) {
    facetax = std::exp(0.522e0_rt*state.eta - 1.563_rt);
  }
  faceta  = 1.0_rt + facetax;
  ocompt = 6.65205e-25_rt/(fact * faceta) * state.xne/state.rho;
  orad += ocompt;

  // cutoff radiative opacity when 4kt/hbar is less than the plasma
  // frequency
  Real tcut = con2 * std::sqrt(state.xne);
  if (state.T < tcut) {
    if (tcut > 200.0_rt*state.T) {
      orad = orad * 2.658e86_rt;
    } else {
      cutfac = std::exp(tcut/state.T - 1.0_rt);
      orad = orad * cutfac;
    }
  }

  // fudge molecular opacity for low temps
  Real xkf = t7peek * state.rho * std::pow(state.T * 1.0e-7_rt, 4);
  orad = xkf * orad/(xkf + orad);


  // conductivity section:
  // drel is the dividing line between nondegenerate and degenerate regions,
  // taken from clayton eq. 2-34. if the density is larger than drel, then
  // use the degenerate expressions. if the density is smaller than
  // drelim, use the non-degenerate formulas. in between drel and drelim,
  // apply a smooth blending of the two.

  Real dlog10 = std::log10(state.rho);

  Real drel = 2.4e-7_rt * zbar/abar * state.T * std::sqrt(state.T);
  if (state.T <= 1.0e5_rt) {
    drel = drel * 15.0_rt;
  }
  Real drel10 = std::log10(drel);
  Real drelim = drel10 + 1.0_rt;

  // from iben apj 196 525 1975 for non-degenerate regimes
  if (dlog10 < drelim) {
    Real zdel = state.xne/(avo*t6*std::sqrt(t6));
    Real zdell10 = std::log10(zdel);
    Real eta0 = std::exp(-1.20322_rt + twoth * std::log(zdel));
    eta02 = eta0*eta0;

    // thpl factor
    Real thpl;
    if (zdell10 < 0.645_rt) {
      thpl = -7.5668_rt + std::log(zdel * (1.0_rt + 0.024417_rt*zdel));
    } else {
      if (zdell10 < 2.5_rt) {
        thpl = -7.58110_rt + std::log(zdel*(1.0_rt + 0.02804_rt*zdel));
        if (zdell10 >= 2.0_rt) {
          Real thpla = thpl;
          thpl = -11.0742_rt + std::log(zdel*zdel * (1.0_rt + 9.376_rt/eta02));
          thpl = 2.0_rt*((2.5_rt-zdell10)*thpla + (zdell10-2.0_rt)*thpl);
        }
      } else {
        thpl = -11.0742_rt + std::log(zdel*zdel * (1.0_rt + 9.376_rt/eta02));
      }
    }

    // pefac and walf factors
    Real pefac;
    if (zdell10 < 2.0_rt) {
      pefac = 1.0_rt + 0.021876_rt*zdel;
      if (zdell10 > 1.5_rt) {
        Real pefacal = std::log(pefac);
        Real pefacl = std::log(0.4_rt * eta0 + 1.64496_rt/eta0);
        Real cfac1 = 2.0_rt - zdell10;
        Real cfac2 = zdell10 - 1.5_rt;
        pefac = std::exp(2.0_rt * (cfac1*pefacal + cfac2*pefacl));
      }
    } else {
      pefac = 0.4_rt * eta0 + 1.64496_rt/eta0;
    }

    Real dnefac;
    if (zdel < 40.0_rt) {
      dnefac = 1.0_rt + zdel * (3.4838e-4_rt * zdel - 2.8966e-2_rt);
    } else {
      dnefac = 1.5_rt/eta0 * (1.0_rt - 0.8225_rt/eta02);
    }
    wpar2 = 9.24735e-3_rt * zdel *
      (state.rho*avo*(w[3]+w[4]+w[5])/state.xne + dnefac)/(std::sqrt(t6)*pefac);
    walf = 0.5_rt * std::log(wpar2);
    walf10 = 0.5_rt * std::log10(wpar2);

    // thx, thy and thc factors
    Real thx;
    if (walf10 <= -3.0_rt) {
      thx = std::exp(2.413_rt - 0.124_rt*walf);
    } else if (walf10 <= -1.0_rt) {
      thx = std::exp(0.299_rt - walf*(0.745_rt + 0.0456_rt*walf));
    } else {
      thx = std::exp(0.426_rt - 0.558_rt*walf);
    }

    Real thy;
    if (walf10 <= -3.0_rt) {
      thy = std::exp(2.158_rt - 0.111_rt*walf);
    } else if (walf10 <= 0.0_rt) {
      thy = std::exp(0.553_rt - walf*(0.55_rt + 0.0299_rt*walf));
    } else {
      thy = std::exp(0.553_rt - 0.6_rt*walf);
    }

    Real thc;
    if (walf10 <= -2.5_rt) {
      thc = std::exp(2.924_rt - 0.1_rt*walf);
    } else if (walf10 <= 0.5_rt) {
      thc = std::exp(1.6740_rt - walf*(0.511_rt + 0.0338_rt*walf));
    } else {
      thc = std::exp(1.941_rt - 0.785_rt*walf);
    }

    oh = (xh*thx + xhe*thy + w[5]*third*thc) / (t6*std::exp(thpl));
  }

  // from yakovlev & urpin soviet astro 1980 24 303 and
  // potekhin et al. 1997 aa 323 415 for degenerate regimes
  if (dlog10 > drel10) {
    Real xmas = meff * std::pow(state.xne, third);
    Real ymas = std::sqrt(1.0_rt + xmas*xmas);
    Real wfac = weid * state.T/ymas * state.xne;
    Real cint = 1.0_rt;

    // ion-electron collision frequency and the thermal conductivity
    Real vie = iec * zbar * ymas * cint;
    Real cie = wfac/vie;

    // electron-electron collision frequency and thermal conductivity
    Real tpe = xec * std::sqrt(state.xne/ymas);
    Real yg = rt3 * tpe/state.T;
    Real xrel = 1.009_rt * std::pow(zbar/abar * state.rho * 1.0e-6_rt, third);
    Real beta2 = xrel*xrel/(1.0_rt + xrel*xrel);
    Real jy = (1.0_rt + 6.0_rt/(5.0_rt*xrel*xrel) + 2.0_rt/(5.0_rt*xrel*xrel*xrel*xrel))
      * ( yg*yg*yg / (3.0_rt * std::pow(1.0_rt + 0.07414_rt * yg, 3))
          * std::log((2.81_rt - 0.810_rt*beta2 + yg)/yg)
          + std::pow(M_PI, 5/6.0_rt) * std::pow(yg/(13.91_rt + yg), 4));
    Real vee = 0.511_rt * state.T*state.T * xmas/(ymas*ymass) * std::sqrt(xmas/ymas) * jy;
    Real cee = wfac/vee;

    // total electron thermal conductivity and conversion to an opacity
    Real ov1 = cie * cee/(cee + cie);
    ov = k2c/(ov1*state % rho) * state.T*state.T*state.T;
  }

  // blend the opacities in the intermediate region
  Real ocond;
  if (dlog10 <= drel10) {
    ocond = oh;
  } else if (dlog10 > drel10 && dlog10 < drelim) {
    Real x = state.rho;
    Real x1 = std::pow(10.0_rt, drel10);
    Real x2 = std::pow(10.0_rt, drelim);
    Real alfa = (x-x2)/(x1-x2);
    Real beta = (x-x1)/(x2-x1);
    ocond = alfa*oh + beta*ov;
  } else if (dlog10 >= drelim) {
    ocond = ov;
  }

  // total opacity
  opac = orad * ocond / (ocond + orad);

  state % conductivity = k2c * state.T*state.T*state.T / (opac * state.rho);
}
#endif
