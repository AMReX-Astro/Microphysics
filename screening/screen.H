class screen_factors_t {

public:

  amrex::Real z1 = -1;
  amrex::Real z2 = -1;
  amrex::Real a1 = -1;
  amrex::Real a2 = -1;

  // zs13    = (z1+z2)**(1./3.)
  // zhat    = combination of z1 and z2 raised to the 5/3 power
  // zhat2   = combination of z1 and z2 raised to the 5/12 power
  // lzav    = log of effective charge
  // aznut   = combination of a1,z1,a2,z2 raised to 1/3 power

  amrex::Real zs13 = 0.0;
  amrex::Real zs13inv = 0.0;
  amrex::Real zhat = 0.0;
  amrex::Real zhat2 = 0.0;
  amrex::Real lzav = 0.0;
  amrex::Real aznut = 0.0;

  bool validate_nuclei(const Real z1_pass, const Real z2_pass,
                       const Real a1_pass, const Real a2_pass) {
    // a simple function for unit testing / debug runs to
    // ensure that we are accessing the proper screening info

    return (z1_pass == z1) &&
           (z2_pass == z2) &&
           (a1_pass == a1) &&
           (a2_pass == a2);
  }
};


struct plasma_state_t {

  amrex::Real qlam0z;
  amrex::Real qlam0zdt;
  //amrex::Real qlam0zdd;

  amrex::Real taufac;
  amrex::Real taufacdt;

  amrex::Real aa;
  amrex::Real daadt;
  //amrex::Real daadd;
};


namespace scrn {
  const amrex::Real fact = 1.25992104989487e0_rt;
  const amrex::Real co2        = (1.0_rt/3.0_rt) * 4.248719e3_rt;
  const amrex::Real gamefx     = 0.3e0_rt;
  const amrex::Real gamefs     = 0.8e0_rt;
  const amrex::Real h12_max    = 300.e0_rt;
};


AMREX_FORCE_INLINE
void
fill_plasma_state(plasma_state_t state, const Real temp, const Real dens, Real* y) {

  Real sum = 0.0_rt;
  for (int n = 0; n < NumSpec; n++) {
    sum += y[n];
  }
  Real abar = 1.0_rt / sum;
  Real ytot = sum;

  sum = 0.0_rt;
  Real sum2 = 0.0_rt;
  for (int n = 0; n < NumSpec; n++) {
    sum += zion[n]*y[n];
    sum2 += zion[n]*zion[n]*y[n];
  }

  Real zbar = sum * abar;
  Real z2bar = sum2 * abar;

  Real rr = dens * ytot;
  Real tempi = 1.0_rt / temp;
  Real dtempi = -tempi * tempi;
  Real deni = 1.0_rt / dens;

  Real pp = std::sqrt(rr*tempi*(z2bar + zbar));
  Real qq = 0.5_rt/pp *(z2bar + zbar);
  Real dppdt = qq*rr*dtempi;
  //Real dppdd = qq * ytot * tempi;

  state.qlam0z = 1.88e8_rt * tempi * pp;
  state.qlam0zdt = 1.88e8_rt * (dtempi*pp + tempi*dppdt);
  //state.qlam0zdd = 1.88e8_rt * tempi * dppdd;

  state.taufac = co2 * std::pow(tempi, 1.0_rt/3.0_rt);
  state.taufacdt = -(1.0_rt/3.0_rt) * state % taufac * tempi;

  qq = rr * zbar;
  xni = std::pow(qq, 1.0_rt/3.0_rt);
  //dxnidd = (1.0_rt/3.0_rt) * xni * deni;

  state.aa = 2.27493e5_rt * tempi * xni;
  state.daadt = 2.27493e5_rt * dtempi * xni;
  //state.daadd = 2.27493e5_rt * tempi * dxnidd;
}
