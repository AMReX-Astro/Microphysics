#include <actual_network.H>

#ifndef _eos_H_
#define _eos_H_

using namespace amrex;

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



#endif
