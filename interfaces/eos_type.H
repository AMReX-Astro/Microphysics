#ifndef EOS_TYPE_H
#define EOS_TYPE_H

#include <AMReX.H>
#include <network.H>

struct eos_base_t {};

struct eos_t:eos_base_t {
    amrex::Real rho{};
    amrex::Real T{};
    amrex::Real p{};
    amrex::Real e{};
    amrex::Real h{};
    amrex::Real s{};
    amrex::Real xn[NumSpec]{};
#if NAUX_NET > 0
    amrex::Real aux[NumAux]{};
#endif
    amrex::Real dpdT{};
    amrex::Real dpdr{};
    amrex::Real dedT{};
    amrex::Real dedr{};
    amrex::Real dhdT{};
    amrex::Real dhdr{};
    amrex::Real dsdT{};
    amrex::Real dsdr{};
    amrex::Real dpde{};
    amrex::Real dpdr_e{};

    amrex::Real cv{};
    amrex::Real cp{};
    amrex::Real xne{};
    amrex::Real xnp{};
    amrex::Real eta{};
    amrex::Real pele{};
    amrex::Real ppos{};
    amrex::Real mu{};
    amrex::Real mu_e{};
    amrex::Real y_e{};
    amrex::Real gam1{};
    amrex::Real cs{};

    amrex::Real abar{};
    amrex::Real zbar{};

#ifdef EXTRA_THERMO
    amrex::Real dpdA{};
    amrex::Real dpdZ{};
    amrex::Real dedA{};
    amrex::Real dedZ{};
#endif
    amrex::Real conductivity{};

#ifdef NSE_NET
    amrex::Real mu_p{};   //chemical potential of proton
    amrex::Real mu_n{};   //chemical potential of neutron
#endif
};

struct eos_re_t:eos_base_t {
    amrex::Real rho{};
    amrex::Real T{};
    amrex::Real e{};
    amrex::Real xn[NumSpec]{};
#if NAUX_NET > 0
    amrex::Real aux[NumAux]{};
#endif
    amrex::Real dedT{};
    amrex::Real dedr{};

    amrex::Real cv{};

    amrex::Real mu{};
    amrex::Real mu_e{};
    amrex::Real y_e{};
    amrex::Real eta{};

    amrex::Real abar{};
    amrex::Real zbar{};

#ifdef EXTRA_THERMO
    amrex::Real dedA{};
    amrex::Real dedZ{};
#endif

#ifdef NSE_NET
    amrex::Real mu_p{};   //chemical potential of proton
    amrex::Real mu_n{};   //chemical potential of neutron
#endif
};

struct eos_rep_t:eos_base_t {
    amrex::Real rho{};
    amrex::Real T{};
    amrex::Real e{};
    amrex::Real p{};
    amrex::Real xn[NumSpec]{};
#if NAUX_NET > 0
    amrex::Real aux[NumAux]{};
#endif
    amrex::Real dedT{};
    amrex::Real dedr{};

    amrex::Real dpdT{};
    amrex::Real dpdr{};
    amrex::Real dpdr_e{};
    amrex::Real dpde{};

    amrex::Real cv{};
    amrex::Real cp{};
    amrex::Real gam1{};
    amrex::Real cs{};

    amrex::Real mu{};
    amrex::Real mu_e{};
    amrex::Real y_e{};
    amrex::Real eta{};

    amrex::Real abar{};
    amrex::Real zbar{};

#ifdef EXTRA_THERMO
    amrex::Real dpdA{};
    amrex::Real dpdZ{};
    amrex::Real dedA{};
    amrex::Real dedZ{};
#endif
  
#ifdef NSE_NET
    amrex::Real mu_p{};   //chemical potential of proton
    amrex::Real mu_n{};   //chemical potential of neutron
#endif
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

template <typename T, typename Enable = void>
struct has_conductivity
    : std::false_type {};

template <typename T>
struct has_conductivity<T, typename std::enable_if<(sizeof(((T*)0)->conductivity) > 0)>::type>
    : std::true_type {};

template <typename T, typename Enable = void>
struct has_base_variables
    : std::false_type {};

template <typename T>
struct has_base_variables<T, typename std::enable_if<(sizeof(((T*)0)->rho) > 0)>::type>
    : std::true_type {};

template <typename T>
inline
std::ostream& print_state (std::ostream& o, T const& eos_state)
{
    if constexpr (has_base_variables<T>::value) {
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
    }
    if constexpr (has_pressure<T>::value) {
        o << "p = " << eos_state.p << std::endl;
    }
    if constexpr (has_energy<T>::value) {
        o << "e = " << eos_state.e << std::endl;
    }
    if constexpr (has_enthalpy<T>::value) {
        o << "h = " << eos_state.h << std::endl;
    }
    if constexpr (has_entropy<T>::value) {
        o << "s = " << eos_state.s << std::endl;
    }

    if constexpr (has_pressure<T>::value) {
        o << "dpdT = " << eos_state.dpdT << std::endl;
        o << "dpdr = " << eos_state.dpdr << std::endl;
    }
    if constexpr (has_energy<T>::value) {
        o << "dedT = " << eos_state.dedT << std::endl;
        o << "dedr = " << eos_state.dedr << std::endl;
    }
    if constexpr (has_enthalpy<T>::value) {
        o << "dhdT = " << eos_state.dhdT << std::endl;
        o << "dhdr = " << eos_state.dhdr << std::endl;
    }
    if constexpr (has_entropy<T>::value) {
        o << "dsdT = " << eos_state.dsdT << std::endl;
        o << "dsdr = " << eos_state.dsdr << std::endl;
    }
    if constexpr (has_pressure<T>::value) {
        o << "dpde = " << eos_state.dpde << std::endl;
        o << "dpdr_e = " << eos_state.dpdr_e  << std::endl;
    }

    if constexpr (has_base_variables<T>::value) {
        o << "cv = " << eos_state.cv << std::endl;
    }
    if constexpr (has_pressure<T>::value) {
        o << "cp = " << eos_state.cp << std::endl;
    }
    if constexpr (has_xne_xnp<T>::value) {
        o << "xne = " << eos_state.xne << std::endl;
        o << "xnp = " << eos_state.xnp << std::endl;
    }
    if constexpr (has_eta<T>::value) {
        o << "eta = " << eos_state.eta << std::endl;
    }
    if constexpr (has_pele_ppos<T>::value) {
        o << "pele = " << eos_state.pele << std::endl;
        o << "ppos = " << eos_state.ppos << std::endl;
    }
    if constexpr (has_base_variables<T>::value) {
        o << "mu = " << eos_state.mu << std::endl;
        o << "mu_e = " << eos_state.mu_e << std::endl;
        o << "y_e = " << eos_state.y_e << std::endl;
    }
    if constexpr (has_pressure<T>::value) {
        o << "gam1 = " << eos_state.gam1 << std::endl;
        o << "cs = " << eos_state.cs << std::endl;
    }

    if constexpr (has_base_variables<T>::value) {
        o << "abar = " << eos_state.abar << std::endl;
        o << "zbar = " << eos_state.zbar << std::endl;
    }

#ifdef EXTRA_THERMO
    if constexpr (has_pressure<T>::value) {
        o << "dpdA = " << eos_state.dpdA << std::endl;
        o << "dpdZ = " << eos_state.dpdZ << std::endl;
    }
    if constexpr (has_energy<T>::value) {
        o << "dedA = " << eos_state.dedA << std::endl;
        o << "dedZ = " << eos_state.dedZ << std::endl;
    }
#endif
#ifdef CONDUCTIVITY
    if constexpr (has_conductivity<T>::value) {
        o << "conductivity = " << eos_state.conductivity << std::endl;
    }
#endif

    return o;
}

template <typename T, typename U = typename std::enable_if<std::is_same<T, eos_t>::value ||
                                                           std::is_same<T, eos_re_t>::value ||
                                                           std::is_same<T, eos_rep_t>::value
                                                           >::type>
inline
std::ostream& operator<< (std::ostream& o, T const& eos_state)
{
    return print_state(o, eos_state);
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

namespace EOS {

    enum eos_iter {itemp = 1,
                   idens = 2,
                   iener = 3,
                   ienth = 4,
                   ientr = 5,
                   ipres = 6};

}

#endif
