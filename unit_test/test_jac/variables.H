#ifndef _variables_H_
#define _variables_H_

// A module to provide integer indices into the various storage arrays
// for accessing the different variables by name.

#include <AMReX_Vector.H>
#include <network.H>
#include <burn_type.H>

class plot_t {

public:

    int irho = -1;
    int itemp = 1;
    int ispec_old = -1;
    int ijac = -1;

    int n_plot_comps = 0;

    int next_index(const int num) {
        int next = n_plot_comps;
        n_plot_comps += num;
        return next;
    }

};

plot_t init_variables();

void get_varnames(const plot_t p, amrex::Vector<std::string>& names);

#endif
