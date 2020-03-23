#ifndef _constants_H_
#define _constants_H_

namespace C {
  constexpr amrex::Real avo = 6.0221417930e23_rt;
  constexpr amrex::Real c_light = 2.99792458e10_rt;

  constexpr amrex::Real ev2erg  = 1.60217648740e-12_rt;
  constexpr amrex::Real mev2erg = ev2erg*1.0e6_rt;
  constexpr amrex::Real mev2gr  = mev2erg/c_light**2;

  constexpr amrex::Real mn = 1.67492721184e-24_rt;
  constexpr amrex::Real mp = 1.67262163783e-24_rt;
  constexpr amrex::Real me = 9.1093821545e-28_rt;

  // Conversion factor for the nuclear energy generation rate.
  constexpr amrex::Real enuc_conv2 = -avo*c_light*c_light;
}

#endif
