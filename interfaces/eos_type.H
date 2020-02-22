#ifndef _eos_type_H_
#define _eos_type_H_

#include <AMReX.H>
#include <network.H>

struct eos_t {
    amrex::Real rho;
    amrex::Real T;
    amrex::Real p;
    amrex::Real e;
    amrex::Real h;
    amrex::Real s;
    amrex::Real xn[NumSpec];
    amrex::Real aux[NumAux];

    amrex::Real dpdT;
    amrex::Real dpdr;
    amrex::Real dedT;
    amrex::Real dedr;
    amrex::Real dhdT;
    amrex::Real dhdr;
    amrex::Real dsdT;
    amrex::Real dsdr;
    amrex::Real dpde;
    amrex::Real dpdr_e;

    amrex::Real cv;
    amrex::Real cp;
    amrex::Real xne;
    amrex::Real xnp;
    amrex::Real eta;
    amrex::Real pele;
    amrex::Real ppos;
    amrex::Real mu;
    amrex::Real mu_e;
    amrex::Real y_e;
    amrex::Real gam1;
    amrex::Real cs;

    amrex::Real abar;
    amrex::Real zbar;

#ifdef EXTRA_THERMO
    amrex::Real dpdA;
    amrex::Real dpdZ;
    amrex::Real dedA;
    amrex::Real dedZ;
#endif

    amrex::Real conductivity;

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
  extern AMREX_GPU_MANAGED amrex::Real mintemp;
  extern AMREX_GPU_MANAGED amrex::Real maxtemp;
  extern AMREX_GPU_MANAGED amrex::Real mindens;
  extern AMREX_GPU_MANAGED amrex::Real maxdens;
  extern AMREX_GPU_MANAGED amrex::Real minx;
  extern AMREX_GPU_MANAGED amrex::Real maxx;
  extern AMREX_GPU_MANAGED amrex::Real minye;
  extern AMREX_GPU_MANAGED amrex::Real maxye;
  extern AMREX_GPU_MANAGED amrex::Real mine;
  extern AMREX_GPU_MANAGED amrex::Real maxe;
  extern AMREX_GPU_MANAGED amrex::Real minp;
  extern AMREX_GPU_MANAGED amrex::Real maxp;
  extern AMREX_GPU_MANAGED amrex::Real mins;
  extern AMREX_GPU_MANAGED amrex::Real maxs;
  extern AMREX_GPU_MANAGED amrex::Real minh;
  extern AMREX_GPU_MANAGED amrex::Real maxh;
}

namespace EOSData
{
  bool initialized;
  AMREX_GPU_MANAGED amrex::Real mintemp;
  AMREX_GPU_MANAGED amrex::Real maxtemp;
  AMREX_GPU_MANAGED amrex::Real mindens;
  AMREX_GPU_MANAGED amrex::Real maxdens;
  AMREX_GPU_MANAGED amrex::Real minx;
  AMREX_GPU_MANAGED amrex::Real maxx;
  AMREX_GPU_MANAGED amrex::Real minye;
  AMREX_GPU_MANAGED amrex::Real maxye;
  AMREX_GPU_MANAGED amrex::Real mine;
  AMREX_GPU_MANAGED amrex::Real maxe;
  AMREX_GPU_MANAGED amrex::Real minp;
  AMREX_GPU_MANAGED amrex::Real maxp;
  AMREX_GPU_MANAGED amrex::Real mins;
  AMREX_GPU_MANAGED amrex::Real maxs;
  AMREX_GPU_MANAGED amrex::Real minh;
  AMREX_GPU_MANAGED amrex::Real maxh;
}

#endif
