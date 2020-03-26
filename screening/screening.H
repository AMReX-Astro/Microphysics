#include <screen.H>

#ifndef _screening_H_
#define _screening_H_

void screening_init();

void screening_finalize();

void add_screening_factor(const int i,
                          const amrex::Real z1, const amrex::Real a1,
                          const amrex::Real z2, const amrex::Real a2);

AMREX_GPU_HOST_DEVICE
void screen5(const plasma_state_t state,
             const int jscreen,
             amrex::Real& scor, amrex::Real& scordt, amrex::Real& scordd);


#endif
