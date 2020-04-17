#ifndef _actual_eos_data_H_
#define _actual_eos_data_H_

#include <AMReX.H>
#include <AMReX_REAL.H>

extern AMREX_GPU_MANAGED amrex::Real gamma_const;
extern AMREX_GPU_MANAGED amrex::Real K_const;
extern AMREX_GPU_MANAGED amrex::Real mu_e;
extern AMREX_GPU_MANAGED int  polytrope;
extern AMREX_GPU_MANAGED amrex::Real gm1;
extern AMREX_GPU_MANAGED amrex::Real polytrope_index;

#endif
