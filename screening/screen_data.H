#ifndef _screen_data_H_
#define _screen_data_H_

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_REAL.H>
#include <network_properties.H>

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

      amrex::Real zs13 = 0.0;
      amrex::Real zs13inv = 0.0;
      amrex::Real zhat = 0.0;
      amrex::Real zhat2 = 0.0;
      amrex::Real lzav = 0.0;
      amrex::Real aznut = 0.0;

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

    extern AMREX_GPU_MANAGED amrex::GpuArray<screen_factors_t, NSCREEN> scn_facs;

};

#endif
