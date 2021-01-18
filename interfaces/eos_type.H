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
#if NAUX_NET > 0
    amrex::Real aux[NumAux];
#endif
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

template <typename T, typename Enable = void>
struct has_energy
    : std::false_type {};

template <typename T>
struct has_energy<T, typename std::enable_if<(sizeof(((T*)0)->e) > 0)>::type>
    : std::true_type {};

template <typename T, typename Enable = void>
struct has_enthalpy
    : std::false_type {};

template <typename T>
struct has_enthalpy<T, typename std::enable_if<(sizeof(((T*)0)->h) > 0)>::type>
    : std::true_type {};

template <typename T, typename Enable = void>
struct has_entropy
    : std::false_type {};

template <typename T>
struct has_entropy<T, typename std::enable_if<(sizeof(((T*)0)->s) > 0)>::type>
    : std::true_type {};

template <typename T, typename Enable = void>
struct has_pressure
    : std::false_type {};

template <typename T>
struct has_pressure<T, typename std::enable_if<(sizeof(((T*)0)->p) > 0)>::type>
    : std::true_type {};

template <typename T, typename Enable = void>
struct has_pele_ppos
    : std::false_type {};

template <typename T>
struct has_pele_ppos<T, typename std::enable_if<(sizeof(((T*)0)->pele) > 0)>::type>
    : std::true_type {};

template <typename T, typename Enable = void>
struct has_xne_xnp
    : std::false_type {};

template <typename T>
struct has_xne_xnp<T, typename std::enable_if<(sizeof(((T*)0)->xne) > 0)>::type>
    : std::true_type {};

template <typename T, typename Enable = void>
struct has_eta
    : std::false_type {};

template <typename T>
struct has_eta<T, typename std::enable_if<(sizeof(((T*)0)->eta) > 0)>::type>
    : std::true_type {};

inline
std::ostream& operator<< (std::ostream& o, eos_t const& eos_state)
{
  o << "rho = " << eos_state.rho << std::endl;
  o << "T =   " << eos_state.T << std::endl;
  o << "xn = ";
  for (int n = 0; n < NumSpec; ++n) {
    o << eos_state.xn[n] << " ";
  }
  o << std::endl;
#if NAUX_NET > 0
  o << "aux = ";
  for (int n = 0; n < NumAux; ++n) {
    o << eos_state.aux[n] << " ";
  }
  o << std::endl;
#endif
  o << "p = " << eos_state.p << std::endl;
  o << "e = " << eos_state.e << std::endl;
  o << "h = " << eos_state.h << std::endl;
  o << "s = " << eos_state.s << std::endl;

  o << "dpdT = " << eos_state.dpdT << std::endl;
  o << "dpdr = " << eos_state.dpdr << std::endl;
  o << "dedT = " << eos_state.dedT << std::endl;
  o << "dedr = " << eos_state.dedr << std::endl;
  o << "dhdT = " << eos_state.dhdT << std::endl;
  o << "dhdr = " << eos_state.dhdr << std::endl;
  o << "dsdT = " << eos_state.dsdT << std::endl;
  o << "dsdr = " << eos_state.dsdr << std::endl;
  o << "dpde = " << eos_state.dpde << std::endl;
  o << "dpdr_e = " << eos_state.dpdr_e  << std::endl;

  o << "cv = " << eos_state.cv << std::endl;
  o << "cp = " << eos_state.cp << std::endl;
  o << "xne = " << eos_state.xne << std::endl;
  o << "xnp = " << eos_state.xnp << std::endl;
  o << "eta = " << eos_state.eta << std::endl;
  o << "pele = " << eos_state.pele << std::endl;
  o << "ppos = " << eos_state.ppos << std::endl;
  o << "mu = " << eos_state.mu << std::endl;
  o << "mu_e = " << eos_state.mu_e << std::endl;
  o << "y_e = " << eos_state.y_e << std::endl;
  o << "gam1 = " << eos_state.gam1 << std::endl;
  o << "cs = " << eos_state.cs << std::endl;

  o << "abar = " << eos_state.abar << std::endl;
  o << "zbar = " << eos_state.zbar << std::endl;

#ifdef EXTRA_THERMO
  o << "dpdA = " << eos_state.dpdA << std::endl;
  o << "dpdZ = " << eos_state.dpdZ << std::endl;
  o << "dedA = " << eos_state.dedA << std::endl;
  o << "dedZ = " << eos_state.dedZ << std::endl;
#endif

  o << "conductivity = " << eos_state.conductivity << std::endl;

  return o;
}

enum eos_input_t {eos_input_rt = 0,
                  eos_input_rh,
                  eos_input_tp,
                  eos_input_rp,
                  eos_input_re,
                  eos_input_ps,
                  eos_input_ph,
                  eos_input_th};

// these are used to allow for a generic interface to the
// root finding

enum eos_iter {itemp = 1,
               idens = 2,
               iener = 3,
               ienth = 4,
               ientr = 5,
               ipres = 6};

#endif
