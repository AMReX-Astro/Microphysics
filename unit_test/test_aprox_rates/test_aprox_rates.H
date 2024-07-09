#ifndef TEST_APROX_RATES_H
#define TEST_APROX_RATES_H

#include <extern_parameters.H>
#include <variables.H>

void main_main();

void aprox_rates_test(const Box& bx,
                      const amrex::Real dlogrho, const amrex::Real dlogT, const amrex::Real dmetal,
                      const plot_t& vars,
                      amrex::Array4<amrex::Real> const sp);

void aprox_rates_extra_c12ag(const Box& bx,
                             const amrex::Real dlogrho, const amrex::Real dlogT, const amrex::Real dmetal,
                             const plot_t& vars,
                             amrex::Array4<amrex::Real> const sp);

#endif
