#ifndef _eos_type_H_
#define _eos_type_H_

#include <AMReX.H>
#include <network.H>

struct eos_t {
    Real rho;
    Real T;
    Real p;
    Real e;
    Real h;
    Real s;
    Real xn[NumSpec];
    Real aux[NumAux];

    Real dpdT;
    Real dpdr;
    Real dedT;
    Real dedr;
    Real dhdT;
    Real dhdr;
    Real dsdT;
    Real dsdr;
    Real dpde;
    Real dpdr_e;

    Real cv;
    Real cp;
    Real xne;
    Real xnp;
    Real eta;
    Real pele;
    Real ppos;
    Real mu;
    Real mu_e;
    Real y_e;
    Real gam1;
    Real cs;

    Real abar;
    Real zbar;

#ifdef EXTRA_THERMO
    Real dpdA;
    Real dpdZ;
    Real dedA;
    Real dedZ;
#endif

    Real conductivity;

};

enum eos_input_t {eos_input_rt = 0,
                  eos_input_rh,
                  eos_input_tp,
                  eos_input_rp,
                  eos_input_re,
                  eos_input_ps,
                  eos_input_ph,
                  eos_input_th};

// store some global limits set at initialization time

namespace EOSData
{
  extern bool initialized;
  extern AMREX_GPU_MANAGED Real mintemp;
  extern AMREX_GPU_MANAGED Real maxtemp;
  extern AMREX_GPU_MANAGED Real mindens;
  extern AMREX_GPU_MANAGED Real maxdens;
  extern AMREX_GPU_MANAGED Real minx;
  extern AMREX_GPU_MANAGED Real maxx;
  extern AMREX_GPU_MANAGED Real minye;
  extern AMREX_GPU_MANAGED Real maxye;
  extern AMREX_GPU_MANAGED Real mine;
  extern AMREX_GPU_MANAGED Real maxe;
  extern AMREX_GPU_MANAGED Real minp;
  extern AMREX_GPU_MANAGED Real maxp;
  extern AMREX_GPU_MANAGED Real mins;
  extern AMREX_GPU_MANAGED Real maxs;
  extern AMREX_GPU_MANAGED Real minh;
  extern AMREX_GPU_MANAGED Real maxh;
}

namespace EOSData
{
  bool initialized;
  AMREX_GPU_MANAGED Real mintemp;
  AMREX_GPU_MANAGED Real maxtemp;
  AMREX_GPU_MANAGED Real mindens;
  AMREX_GPU_MANAGED Real maxdens;
  AMREX_GPU_MANAGED Real minx;
  AMREX_GPU_MANAGED Real maxx;
  AMREX_GPU_MANAGED Real minye;
  AMREX_GPU_MANAGED Real maxye;
  AMREX_GPU_MANAGED Real mine;
  AMREX_GPU_MANAGED Real maxe;
  AMREX_GPU_MANAGED Real minp;
  AMREX_GPU_MANAGED Real maxp;
  AMREX_GPU_MANAGED Real mins;
  AMREX_GPU_MANAGED Real maxs;
  AMREX_GPU_MANAGED Real minh;
  AMREX_GPU_MANAGED Real maxh;
}

#endif
