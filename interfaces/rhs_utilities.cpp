#include <rhs_utilities.H>

AMREX_GPU_MANAGED Array3D<Real, 1, Rates::NumRatesFR, 1, 2, 1, RHS::nrattab> RHS::rattab;
AMREX_GPU_MANAGED Array3D<Real, 1, Rates::NumRatesFR, 1, 2, 1, RHS::nrattab> RHS::drattabdt;
AMREX_GPU_MANAGED Array1D<Real, 1, RHS::nrattab> RHS::ttab;
