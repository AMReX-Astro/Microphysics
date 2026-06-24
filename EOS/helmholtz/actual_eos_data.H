#ifndef actual_eos_data_H
#define actual_eos_data_H

#include <AMReX.H>
#include <AMReX_REAL.H>

#if defined(USE_NEW_HELMHOLTZ_TABLE)
#include <fundamental_constants.H>
#endif

namespace helmholtz
{

    using namespace amrex::literals;

    extern AMREX_GPU_MANAGED int do_coulomb;
    extern AMREX_GPU_MANAGED int input_is_constant;

    // for the tables

#if defined(USE_NEW_HELMHOLTZ_TABLE)
    constexpr int imax = 421;
    constexpr int jmax = 161;

    // these are log_2 using the fastlog approximation
    constexpr amrex::Real tlo   = 9.953125_rt;
    constexpr amrex::Real thi   = 36.4551915228366851806640625_rt;
    constexpr amrex::Real tstp  = (thi - tlo) / (static_cast<amrex::Real>(jmax-1));
    constexpr amrex::Real tstpi = 1.0e0_rt / tstp;
    constexpr amrex::Real dlo   = -33.2820130816_rt;
    constexpr amrex::Real dhi   = 36.4551915228366851806640625_rt;
    constexpr amrex::Real dstp  = (dhi - dlo) / (static_cast<amrex::Real>(imax-1));
    constexpr amrex::Real dstpi = 1.0e0_rt / dstp;

#else

    constexpr int imax = 541;
    constexpr int jmax = 201;

    constexpr amrex::Real tlo   = 3.0e0_rt;
    constexpr amrex::Real thi   = 13.0e0_rt;
    constexpr amrex::Real tstp  = (thi - tlo) / (static_cast<amrex::Real>(jmax-1));
    constexpr amrex::Real tstpi = 1.0e0_rt / tstp;
    constexpr amrex::Real dlo   = -12.0e0_rt;
    constexpr amrex::Real dhi   = 15.0e0_rt;
    constexpr amrex::Real dstp  = (dhi - dlo) / (static_cast<amrex::Real>(imax-1));
    constexpr amrex::Real dstpi = 1.0e0_rt / dstp;
#endif

    extern AMREX_GPU_MANAGED amrex::Real d[imax];
    extern AMREX_GPU_MANAGED amrex::Real t[jmax];

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

#if defined(USE_NEW_HELMHOLTZ_TABLE)
    constexpr amrex::Real h = C::hplanck;
    constexpr amrex::Real avo_eos = C::n_A;
    constexpr amrex::Real kerg = C::k_B;
    constexpr amrex::Real amu = C::m_u;
#else
    // 2006 CODATA physical constants
    constexpr amrex::Real h = 6.6260689633e-27;
    constexpr amrex::Real avo_eos = 6.0221417930e23;
    constexpr amrex::Real kerg = 1.380650424e-16;
    constexpr amrex::Real amu = 1.66053878283e-24;
#endif

}

#endif
