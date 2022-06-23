#ifndef screen_data_H
#define screen_data_H

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <network.H>
#include <microphysics_math.H>

namespace scrn {
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
      // ztilde  = effective ion radius factor for a MCP

      amrex::Real zs13 = 0.0;
      amrex::Real zs13inv = 0.0;
      amrex::Real zhat = 0.0;
      amrex::Real zhat2 = 0.0;
      amrex::Real lzav = 0.0;
      amrex::Real aznut = 0.0;
      amrex::Real ztilde = 0.0;

      bool validate_nuclei(const amrex::Real z1_pass, const amrex::Real a1_pass,
                           const amrex::Real z2_pass, const amrex::Real a2_pass) {
        // a simple function for unit testing / debug runs to
        // ensure that we are accessing the proper screening info

        return (z1_pass == z1) &&
               (z2_pass == z2) &&
               (a1_pass == a1) &&
               (a2_pass == a2);
      }
    };

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    constexpr
    screen_factors_t calculate_screen_factor (amrex::Real z1, amrex::Real a1,
                                              amrex::Real z2, amrex::Real a2)
    {
        screen_factors_t scn_fac;

        scn_fac.z1 = z1;
        scn_fac.a1 = a1;
        scn_fac.z2 = z2;
        scn_fac.a2 = a2;

        scn_fac.zs13 = gcem::pow(z1 + z2, 1.0_rt / 3.0_rt);

        scn_fac.zs13inv = 1.0_rt / scn_fac.zs13;

        scn_fac.zhat = gcem::pow(z1 + z2, 5.0_rt / 3.0_rt) -
                       gcem::pow(z1, 5.0_rt/3.0_rt) -
                       gcem::pow(z2, 5.0_rt/3.0_rt);

        scn_fac.zhat2 = gcem::pow(z1 + z2, 5.0_rt / 12.0_rt) -
                        gcem::pow(z1, 5.0_rt/12.0_rt) -
                        gcem::pow(z2, 5.0_rt/12.0_rt);

        scn_fac.lzav = (5.0_rt / 3.0_rt) * gcem::log(z1 * z2 / (z1 + z2));

        scn_fac.aznut = gcem::pow(z1 * z1 * z2 * z2 * a1 * a2 / (a1 + a2),
                                  1.0_rt / 3.0_rt);

        scn_fac.ztilde = 0.5_rt * (gcem::pow(z1, 1.0_rt / 3.0_rt) +
                                   gcem::pow(z2, 1.0_rt / 3.0_rt));

        return scn_fac;
    }

#if NUMSCREEN > 0
    extern AMREX_GPU_MANAGED amrex::GpuArray<screen_factors_t, NSCREEN> scn_facs;
#endif
}

#endif
