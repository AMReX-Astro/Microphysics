#include <aprox_rates_data.H>

using namespace amrex;

AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,6> rv = {6.0_rt, 7.0_rt, 8.0_rt, 9.0_rt, 10.0_rt, 11.0_rt};
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,14> tv = {1.0_rt,2.0_rt,3.0_rt,4.0_rt,5.0_rt,6.0_rt,7.0_rt,8.0_rt,9.0_rt,10.0_rt,11.0_rt,12.0_rt,13.0_rt,14.0_rt};
AMREX_GPU_MANAGED amrex::Array3D<amrex::Real,1,2,1,6,1,14> datn;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,4> rfdm;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,4> rfd0;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,4> rfd1;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,4> rfd2;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,12> tfdm;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,12> tfd0;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,12> tfd1;
AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,12> tfd2;

