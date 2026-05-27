#ifndef actual_eos_data_H
#define actual_eos_data_H

#include <AMReX.H>
#include <AMReX_REAL.H>

inline AMREX_GPU_MANAGED amrex::Real gamma_const;
inline AMREX_GPU_MANAGED amrex::Real K_const;
inline AMREX_GPU_MANAGED amrex::Real mu_e;
inline AMREX_GPU_MANAGED int  polytrope;
inline AMREX_GPU_MANAGED amrex::Real gm1;
inline AMREX_GPU_MANAGED amrex::Real polytrope_index;

#endif
