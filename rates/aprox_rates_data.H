#ifndef _aprox_rates_data_H_
#define _aaprox_rates_data_H_

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_REAL.H>

extern amrex::GpuArray<amrex::Real,6> rv;
extern amrex::GpuArray<amrex::Real,14> tv;
extern amrex::Vector<amrex::Vector<amrex::GpuArray<amrex::Real,6> > > datn;
extern amrex::GpuArray<amrex::Real,4> rfdm;
extern amrex::GpuArray<amrex::Real,4> rfd0;
extern amrex::GpuArray<amrex::Real,4> rfd1;
extern amrex::GpuArray<amrex::Real,4> rfd2;
extern amrex::GpuArray<amrex::Real,12> tfdm;
extern amrex::GpuArray<amrex::Real,12> tfd0;
extern amrex::GpuArray<amrex::Real,12> tfd1;
extern amrex::GpuArray<amrex::Real,12> tfd2;        

#endif
