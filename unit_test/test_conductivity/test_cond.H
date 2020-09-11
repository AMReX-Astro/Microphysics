#ifndef TEST_COND_H
#define TEST_COND_H

#include <extern_parameters.H>
#include <variables.H>

void main_main();

void cond_test_C(const Box& bx,
                 const Real dlogrho, const Real dlogT, const Real dmetal,
                 const plot_t vars,
                 Array4<Real> const sp);

#endif
