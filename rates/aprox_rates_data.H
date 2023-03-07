#ifndef APROX_RATES_DATA_H
#define APROX_RATES_DATA_H

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_REAL.H>

extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,6> rv;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,14> tv;
extern AMREX_GPU_MANAGED amrex::Array3D<amrex::Real,1,2,1,6,1,14> datn;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,4> rfdm;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,4> rfd0;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,4> rfd1;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,4> rfd2;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,12> tfdm;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,12> tfd0;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,12> tfd1;
extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real,1,12> tfd2;        

#endif
