#include <actual_eos_data.H>

AMREX_GPU_MANAGED int helmholtz::do_coulomb;
AMREX_GPU_MANAGED int helmholtz::input_is_constant;

AMREX_GPU_MANAGED amrex::Real helmholtz::d[imax];
AMREX_GPU_MANAGED amrex::Real helmholtz::t[jmax];

AMREX_GPU_MANAGED amrex::Real helmholtz::ttol;
AMREX_GPU_MANAGED amrex::Real helmholtz::dtol;

// for the helmholtz free energy tables
AMREX_GPU_MANAGED amrex::Real helmholtz::f[jmax][imax][9];

// for the pressure derivative with density tables
AMREX_GPU_MANAGED amrex::Real helmholtz::dpdf[jmax][imax][4];

// for chemical potential tables
AMREX_GPU_MANAGED amrex::Real helmholtz::ef[jmax][imax][4];

// for the number density tables
AMREX_GPU_MANAGED amrex::Real helmholtz::xf[jmax][imax][4];

// for storing the differences
AMREX_GPU_MANAGED amrex::Real helmholtz::dt_sav[jmax];
AMREX_GPU_MANAGED amrex::Real helmholtz::dt2_sav[jmax];
AMREX_GPU_MANAGED amrex::Real helmholtz::dti_sav[jmax];
AMREX_GPU_MANAGED amrex::Real helmholtz::dt2i_sav[jmax];

AMREX_GPU_MANAGED amrex::Real helmholtz::dd_sav[imax];
AMREX_GPU_MANAGED amrex::Real helmholtz::dd2_sav[imax];
AMREX_GPU_MANAGED amrex::Real helmholtz::ddi_sav[imax];
AMREX_GPU_MANAGED amrex::Real helmholtz::dd2i_sav[imax];
