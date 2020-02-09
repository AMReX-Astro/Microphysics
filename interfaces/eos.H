#include <network_properties.H>

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



// store some global limits set at initialization time

namespace EOSData
{
  extern bool initialized;
  extern Real mintemp;
  extern Real maxtemp;
  extern Real mindens;
  extern Real maxdens;
  extern Real minx;
  extern Real maxx;
  extern Real minye;
  extern Real maxye;
  extern Real mine;
  extern Real maxe;
  extern Real minp;
  extern Real maxp;
  extern Real mins;
  extern Real maxs;
  extern Real minh;
  extern Real maxh;
}

void eos_cxx_init();

void eos_cxx_finalize();

void eos_cxx(const eos_input_t input, eos_t& state, bool use_raw_inputs=false);

void reset_inputs(const eos_input_t input, eos_t& state, bool& has_been_reset);

inline void reset_rho(eos_t& state, bool& has_been_reset);

inline void reset_T(eos_t& state, bool& has_been_reset);

void reset_e(eos_t& state, bool& has_been_reset);

void reset_h(eos_t& state, bool& has_been_reset);

void reset_s(eos_t& state, bool& has_been_reset);

void reset_p(eos_t& state, bool& has_been_reset);

void eos_reset(eos_t& state, bool& has_been_reset);

#ifndef AMREX_USE_GPU
void check_inputs(const eos_input_t input, eos_t& state);

void check_rho(eos_t& state);

void check_T(eos_t& state);

void check_e(eos_t& state);

void check_h(eos_t& state);

void check_s(eos_t& state);

void check_p(eos_t& state);
#endif

inline void print_state(eos_t& state) {
  std::cout << "DENS = " << state.rho << std::endl;
  std::cout << "TEMP = " << state.T << std::endl;
  std::cout << "X = ";
  for (int n = 0; n < NumSpec; n ++) {
    std::cout << state.xn[n] << " ";
  }
  std::cout << std::endl;

  std::cout << "Y_E = " << state.y_e;
}

#endif
