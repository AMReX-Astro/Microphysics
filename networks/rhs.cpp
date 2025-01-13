#include <actual_network.H>

#ifdef NEW_NETWORK_IMPLEMENTATION
#include <rhs.H>

AMREX_GPU_MANAGED amrex::Array3D<autodiff::dual, 1, Rates::NumRates, 1, 2, 1, RHS::nrattab> RHS::rattab{};
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, RHS::nrattab> RHS::ttab;

#endif
