#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>

#include <variables.H>
#include <network.H>
#include <eos.H>

#include <cmath>

using namespace amrex;

void eos_test_C(const Box& bx,
                const Real dlogrho, const Real dlogT, const Real dmetal,
                const plot_t vars,
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

    eos_xderivs_t eos_xderivs = composition_derivatives(eos_state);

    eos_t eos_state_reference;
    eos_state_reference = eos_state;

    sp(i, j, k, vars.ih) = eos_state.h;
    sp(i, j, k, vars.ie) = eos_state.e;
    sp(i, j, k, vars.ip) = eos_state.p;
    sp(i, j, k, vars.is) = eos_state.s;

    sp(i, j, k, vars.icv) = eos_state.cv;
    sp(i, j, k, vars.icp) = eos_state.cp;
    sp(i, j, k, vars.ine) = eos_state.xne;
    sp(i, j, k, vars.inp) = eos_state.xnp;
    sp(i, j, k, vars.ieta) = eos_state.eta;
    sp(i, j, k, vars.ipele) = eos_state.pele;
    sp(i, j, k, vars.ippos) = eos_state.ppos;
    sp(i, j, k, vars.imu) = eos_state.mu;
    sp(i, j, k, vars.imue) = eos_state.mu_e;
    sp(i, j, k, vars.idpdt) = eos_state.dpdT;
    sp(i, j, k, vars.idpdr) = eos_state.dpdr;
    sp(i, j, k, vars.idedt) = eos_state.dedT;
    sp(i, j, k, vars.idedr) = eos_state.dedr;
    sp(i, j, k, vars.idhdt) = eos_state.dhdT;
    sp(i, j, k, vars.idhdr) = eos_state.dhdr;
    sp(i, j, k, vars.idsdt) = eos_state.dsdT;
    sp(i, j, k, vars.idsdr) = eos_state.dsdr;
    for (int n = 0; n < NumSpec; n++) {
      sp(i, j, k, vars.idpdx + n) = eos_xderivs.dpdX[n];
      sp(i, j, k, vars.idedx + n) = eos_xderivs.dedX[n];
      sp(i, j, k, vars.idhdx + n) = eos_xderivs.dhdX[n];
    }
    sp(i, j, k, vars.igam1) = eos_state.gam1;
    sp(i, j, k, vars.ics) = eos_state.cs;
    sp(i, j, k, vars.iabar) = eos_state.abar;
    sp(i, j, k, vars.izbar) = eos_state.zbar;
    sp(i, j, k, vars.idpda) = eos_state.dpdA;
    sp(i, j, k, vars.idpdz) = eos_state.dpdZ;
    sp(i, j, k, vars.ideda) = eos_state.dedA;
    sp(i, j, k, vars.idedz) = eos_state.dedZ;
    sp(i, j, k, vars.idpde) = eos_state.dpde;
    sp(i, j, k, vars.idpdre) = eos_state.dpdr_e;


    // call EOS using rho, h

    // reset T to give it some work to do
    eos_state.T = 100.0;

    eos(eos_input_rh, eos_state);

    sp(i, j, k, vars.ierr_T_eos_rh) =
        std::abs(eos_state.T - temp_zone)/temp_zone;

    eos_state = eos_state_reference;


    // call EOS using T, p
    if (eos_name == "gamma_law") {
      sp(i, j, k, vars.ierr_rho_eos_tp) = 0.0;
    } else {
      // reset rho to give it some work to do
      eos_state.rho = 1.0;

      eos(eos_input_tp, eos_state);

      sp(i, j, k, vars.ierr_rho_eos_tp) =
        std::abs(eos_state.rho - dens_zone)/dens_zone;

      eos_state = eos_state_reference;
    }

    // call EOS using r, p

    // reset T to give it some work to do
    eos_state.T = 100.0;

    eos(eos_input_rp, eos_state);

    sp(i, j, k, vars.ierr_T_eos_rp) =
        std::abs(eos_state.T - temp_zone)/temp_zone;

    eos_state = eos_state_reference;


    // call EOS using r, e

    // reset T to give it some work to do
    eos_state.T = 100.0;

    eos(eos_input_re, eos_state);

    sp(i, j, k, vars.ierr_T_eos_re) =
        std::abs(eos_state.T - temp_zone)/temp_zone;

    eos_state = eos_state_reference;


    // call EOS using p, s

    // reset T and rho to give it some work to do
    eos_state.T = 100.0;
    eos_state.rho = 1.0;


    // some EOSes don't have physically valid treatments
    // of entropy throughout the entire rho-T plane
    if (eos_name == "multigamma" || eos_state.s <= 0.0) {
      sp(i, j, k, vars.ierr_T_eos_ps) = 0.0;
      sp(i, j, k, vars.ierr_rho_eos_ps) = 0.0;

    } else {
      eos(eos_input_ps, eos_state);

      // store the thermodynamic state
      sp(i, j, k, vars.ierr_T_eos_ps) =
          std::abs(eos_state.T - temp_zone)/temp_zone;
      sp(i, j, k, vars.ierr_rho_eos_ps) =
          std::abs(eos_state.rho - dens_zone)/dens_zone;

    }

    eos_state = eos_state_reference;


    // call EOS using p, h

    // reset T and rho to give it some work to do
    eos_state.T = 100.0;
    eos_state.rho = 1.0;

    eos(eos_input_ph, eos_state);

    sp(i, j, k, vars.ierr_T_eos_ph) =
        std::abs(eos_state.T - temp_zone)/temp_zone;
    sp(i, j, k, vars.ierr_rho_eos_ph) =
        std::abs(eos_state.rho - dens_zone)/dens_zone;

    eos_state = eos_state_reference;


    // call EOS using T, h
    // this doesn't work for all EOSes (where h doesn't depend on T)
    if (eos_name == "gamma_law_general" || eos_name == "multigamma") {
      sp(i, j, k, vars.ierr_rho_eos_th) = 0.0;

    } else {
      // reset rho to give it some work to do -- for helmeos, h is not
      // monotonic, so we only perturb rho slightly here
      eos_state.rho = 0.9 * eos_state.rho;

      eos(eos_input_th, eos_state);

      sp(i, j, k, vars.ierr_rho_eos_th) =
        std::abs(eos_state.rho - dens_zone)/dens_zone;

      eos_state = eos_state_reference;
    }

  });
}
