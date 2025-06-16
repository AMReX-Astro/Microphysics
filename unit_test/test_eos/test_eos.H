#ifndef TEST_EOS_H
#define TEST_EOS_H

#include <extern_parameters.H>
#include <variables.H>
#include <react_util.H>

void main_main();

void eos_test_C(const amrex::Box& bx,
                const amrex::Real dlogrho, const amrex::Real dlogT, const init_t& comp_data,
                const plot_t& vars,
                amrex::Array4<amrex::Real> const sp);

#endif
