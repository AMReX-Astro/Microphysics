#ifndef _actual_eos_data_H_
#define _actual_eos_data_H_

#include <AMReX.H>
#include <AMReX_REAL.H>

namespace helmholtz
{

    extern AMREX_GPU_MANAGED bool do_coulomb;
    extern AMREX_GPU_MANAGED bool input_is_constant;

    // for the tables

    const int imax = 541;
    const int jmax = 201;

    extern AMREX_GPU_MANAGED int itmax;
    extern AMREX_GPU_MANAGED int jtmax;

    extern AMREX_GPU_MANAGED amrex::Real d[imax];
    extern AMREX_GPU_MANAGED amrex::Real t[jmax];

    extern AMREX_GPU_MANAGED amrex::Real tlo;
    extern AMREX_GPU_MANAGED amrex::Real thi;
    extern AMREX_GPU_MANAGED amrex::Real tstp;
    extern AMREX_GPU_MANAGED amrex::Real tstpi;

    extern AMREX_GPU_MANAGED amrex::Real dlo;
    extern AMREX_GPU_MANAGED amrex::Real dhi;
    extern AMREX_GPU_MANAGED amrex::Real dstp;
    extern AMREX_GPU_MANAGED amrex::Real dstpi;

    extern AMREX_GPU_MANAGED amrex::Real ttol;
    extern AMREX_GPU_MANAGED amrex::Real dtol;

    // for the helmholtz free energy tables
    extern AMREX_GPU_MANAGED amrex::Real f[jmax][imax][9];

    // for the pressure derivative with density tables
    extern AMREX_GPU_MANAGED amrex::Real dpdf[jmax][imax][4];

    // for chemical potential tables
    extern AMREX_GPU_MANAGED amrex::Real ef[jmax][imax][4];

    // for the number density tables
    extern AMREX_GPU_MANAGED amrex::Real xf[jmax][imax][4];

    // for storing the differences
    extern AMREX_GPU_MANAGED amrex::Real dt_sav[jmax];
    extern AMREX_GPU_MANAGED amrex::Real dt2_sav[jmax];
    extern AMREX_GPU_MANAGED amrex::Real dti_sav[jmax];
    extern AMREX_GPU_MANAGED amrex::Real dt2i_sav[jmax];

    extern AMREX_GPU_MANAGED amrex::Real dd_sav[imax];
    extern AMREX_GPU_MANAGED amrex::Real dd2_sav[imax];
    extern AMREX_GPU_MANAGED amrex::Real ddi_sav[imax];
    extern AMREX_GPU_MANAGED amrex::Real dd2i_sav[imax];

    // 2006 CODATA physical constants
    const amrex::Real h = 6.6260689633e-27;
    const amrex::Real avo_eos = 6.0221417930e23;
    const amrex::Real kerg = 1.380650424e-16;
    const amrex::Real amu = 1.66053878283e-24;

}

#endif
