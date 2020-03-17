#include <vector>
#include <string>

#ifndef _variables_H_
#define _variables_H_

#include <AMReX_Vector.H>

class plot_t {

public:


  int irho = -1;
  int itemp = -1;
  int ih = -1;
  int ie = -1;
  int ip = -1;
  int is = -1;
  int iconductivity = -1;
  int ispec = -1;

  int n_plot_comps = 0;

  amrex::Vector<std::string> names;

  int next_index(const int num) {
    int next = n_plot_comps;
    n_plot_comps += num;
    return next;
  }

};

plot_t init_variables();


#endif
