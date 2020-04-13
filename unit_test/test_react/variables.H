#ifndef _variables_H_
#define _variables_H_

// A module to provide integer indices into the various storage arrays
// for accessing the different variables by name.

#include <network.H>

enum plot_t
{
    p_irho = 0,
    p_itemp = 1,
    p_ispec = 2,
    p_ispec_old = p_ispec + NumSpec,
    p_irodot = p_ispec_old + NumSpec,
    p_irho_hnuc = p_irodot + NumSpec,
    n_plot_comps
};

#endif
