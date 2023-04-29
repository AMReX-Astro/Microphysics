#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>

#include <variables.H>
#include <network.H>
#include <eos.H>
#include <conductivity.H>

#include <cmath>

using namespace amrex;

void cond_test_C(const Box& bx,
                 const Real dlogrho, const Real dlogT, const Real dmetal,
                 const plot_t& vars,
                 Array4<Real> const sp) {

  const int ih1 = network_spec_index("hydrogen-1");
  const int ihe4 = network_spec_index("helium-4");

  AMREX_PARALLEL_FOR_3D(bx, i, j, k,
  {

    // set the composition -- approximately solar
    Real metalicity = 0.0 + static_cast<Real> (k) * dmetal;

    eos_t eos_state;

    for (int n = 0; n < NumSpec; n++) {
      eos_state.xn[n] = metalicity/(NumSpec - 2);
    }
    eos_state.xn[ih1] = 0.75 - 0.5*metalicity;
    eos_state.xn[ihe4] = 0.25 - 0.5*metalicity;

    Real temp_zone = std::pow(10.0, std::log10(temp_min) + static_cast<Real>(j)*dlogT);
    eos_state.T = temp_zone;

    Real dens_zone = std::pow(10.0, std::log10(dens_min) + static_cast<Real>(i)*dlogrho);
    eos_state.rho = dens_zone;

    // store default state
    sp(i, j, k, vars.irho) = eos_state.rho;
    sp(i, j, k, vars.itemp) = eos_state.T;
    for (int n = 0; n < NumSpec; n++) {
      sp(i, j, k, vars.ispec+n) = eos_state.xn[n];
    }

    // call the EOS using rho, T
    eos(eos_input_rt, eos_state);

    conductivity(eos_state);

    sp(i, j, k, vars.ih) = eos_state.h;
    sp(i, j, k, vars.ie) = eos_state.e;
    sp(i, j, k, vars.ip) = eos_state.p;
    sp(i, j, k, vars.is) = eos_state.s;

    sp(i, j, k, vars.iconductivity) = eos_state.conductivity;

  });

}
