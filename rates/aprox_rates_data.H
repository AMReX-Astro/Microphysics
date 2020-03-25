#ifndef _aprox_rates_data_H_
#define _aaprox_rates_data_H_

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_REAL.H>

extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real,6> rv;
extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real,14> tv;
extern AMREX_GPU_MANAGED amrex::Array3D<amrex::Real,0,1,0,5,0,13> datn;
extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real,4> rfdm;
extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real,4> rfd0;
extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real,4> rfd1;
extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real,4> rfd2;
extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real,12> tfdm;
extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real,12> tfd0;
extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real,12> tfd1;
extern AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real,12> tfd2;        

#endif
