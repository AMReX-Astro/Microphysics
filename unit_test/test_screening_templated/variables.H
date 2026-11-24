#include <string>

#ifndef VARIABLES_H
#define VARIABLES_H

#include <AMReX_Vector.H>
#include <AMReX_Array.H>
#include <actual_network.H>

#ifndef NEW_NETWORK_IMPLEMENTATION
#error "only templated networks can be used with this unit test"
#endif

class plot_t {

public:

  struct index_t {
    int value = -1;
    int dt = -1;
    int aux_value = -1;
    int aux_dt = -1;
  };

  int irho = -1;
  int itemp = -1;
  int ispec = -1;

  amrex::Array1D<index_t, 1, Rates::NumRates> iscn;

  int n_plot_comps = 0;

  int next_index(const int num) {
    int next = n_plot_comps;
    n_plot_comps += num;
    return next;
  }

};

plot_t init_variables(amrex::Vector<std::string>& names);


#endif
