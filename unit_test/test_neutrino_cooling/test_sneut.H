#ifndef TEST_SCREEN_H
#define TEST_SCREEN_H

#include <extern_parameters.H>
#include <variables.H>

void main_main();

void neut_test_C(const Box& bx,
                 const amrex::Real dlogrho, const amrex::Real dlogT, const amrex::Real dmetal,
                 const plot_t& vars,
                 amrex::Array4<amrex::Real> const sp);

#endif
