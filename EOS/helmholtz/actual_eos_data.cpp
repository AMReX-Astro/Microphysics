#include <actual_eos_data.H>

AMREX_GPU_MANAGED bool do_coulomb;
AMREX_GPU_MANAGED bool input_is_constant;

AMREX_GPU_MANAGED int itmax;
AMREX_GPU_MANAGED int jtmax;

AMREX_GPU_MANAGED amrex::Real d[imax];
AMREX_GPU_MANAGED amrex::Real t[jmax];

AMREX_GPU_MANAGED amrex::Real tlo;
AMREX_GPU_MANAGED amrex::Real thi;
AMREX_GPU_MANAGED amrex::Real tstp;
AMREX_GPU_MANAGED amrex::Real tstpi;

AMREX_GPU_MANAGED amrex::Real dlo;
AMREX_GPU_MANAGED amrex::Real dhi;
AMREX_GPU_MANAGED amrex::Real dstp;
AMREX_GPU_MANAGED amrex::Real dstpi;

AMREX_GPU_MANAGED amrex::Real ttol;
AMREX_GPU_MANAGED amrex::Real dtol;

// for the helmholtz free energy tables
AMREX_GPU_MANAGED amrex::Real f[jmax][imax][9];

// for the pressure derivative with density tables
AMREX_GPU_MANAGED amrex::Real dpdf[jmax][imax][4];

// for chemical potential tables
AMREX_GPU_MANAGED amrex::Real ef[jmax][imax][4];

// for the number density tables
AMREX_GPU_MANAGED amrex::Real xf[jmax][imax][4];

// for storing the differences
AMREX_GPU_MANAGED amrex::Real dt_sav[jmax];
AMREX_GPU_MANAGED amrex::Real dt2_sav[jmax];
AMREX_GPU_MANAGED amrex::Real dti_sav[jmax];
AMREX_GPU_MANAGED amrex::Real dt2i_sav[jmax];

AMREX_GPU_MANAGED amrex::Real dd_sav[imax];
AMREX_GPU_MANAGED amrex::Real dd2_sav[imax];
AMREX_GPU_MANAGED amrex::Real ddi_sav[imax];
AMREX_GPU_MANAGED amrex::Real dd2i_sav[imax];
