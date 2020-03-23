#include <aprox_rates_data.H>

using namespace amrex;

amrex::GpuArray<amrex::Real,6> rv = {6.0_rt, 7.0_rt, 8.0_rt, 9.0_rt, 10.0_rt, 11.0_rt};
amrex::GpuArray<amrex::Real,14> tv = {1.0_rt,2.0_rt,3.0_rt,4.0_rt,5.0_rt,6.0_rt,7.0_rt,8.0_rt,9.0_rt,10.0_rt,11.0_rt,12.0_rt,13.0_rt,14.0_rt};
amrex::Vector<amrex::Vector<amrex::GpuArray<amrex::Real,6> > > datn;
amrex::GpuArray<amrex::Real,4> rfdm;
amrex::GpuArray<amrex::Real,4> rfd0;
amrex::GpuArray<amrex::Real,4> rfd1;
amrex::GpuArray<amrex::Real,4> rfd2;
amrex::GpuArray<amrex::Real,12> tfdm;
amrex::GpuArray<amrex::Real,12> tfd0;
amrex::GpuArray<amrex::Real,12> tfd1;
amrex::GpuArray<amrex::Real,12> tfd2;