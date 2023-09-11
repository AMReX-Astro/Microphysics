#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>

#include <variables.H>
#include <network.H>
#include <eos.H>

#include <sneut5.H>

#include <cmath>

using namespace amrex;
using namespace unit_test_rp;

void neut_test_C(const Box& bx,
                 const Real dlogrho, const Real dlogT, const Real dmetal,
                 const plot_t& vars,
                 Array4<Real> const sp) {

  using namespace Species;

  const int ih1 = network_spec_index("hydrogen-1");
  if (ih1 < 0) amrex::Error("Error: ih1 not found");

  const int ihe4 = network_spec_index("helium-4");
  if (ihe4 < 0) amrex::Error("Error: ihe4 not found");


  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k)
  {

    // set the composition -- approximately solar
    Real metalicity = 0.0 + static_cast<Real> (k) * dmetal;

    Real xn[NumSpec];

    // for now... the screening using 1-based indexing
    Array1D<Real, 1, NumSpec> ymass;

    for (int n = 0; n < NumSpec; n++) {
      xn[n] = metalicity / static_cast<Real>(NumSpec - 2);
    }
    xn[ih1] = 0.75_rt - 0.5_rt * metalicity;
    xn[ihe4] = 0.25_rt - 0.5_rt * metalicity;

    Real temp_zone = std::pow(10.0, std::log10(temp_min) + static_cast<Real>(j)*dlogT);

    Real dens_zone = std::pow(10.0, std::log10(dens_min) + static_cast<Real>(i)*dlogrho);

    // store default state
    sp(i, j, k, vars.irho) = dens_zone;
    sp(i, j, k, vars.itemp) = temp_zone;
    for (int n = 0; n < NumSpec; n++) {
      sp(i, j, k, vars.ispec+n) = xn[n];
    }

    // compute abar and zbar

    Real abar = 0.0;
    for (int n = 0; n < NumSpec; n++) {
      abar += xn[n] / aion[n];
    }
    abar = 1.0_rt / abar;

    Real zbar = 0.0;
    for (int n = 0; n < NumSpec; n++) {
      zbar += zion[n] * xn[n] / aion[n];
    }
    zbar *= abar;

    Real snu;
    Real dsnudt;
    Real dsnudd;
    Real dsnuda;
    Real dsnudz;

    sneut5(temp_zone, dens_zone, abar, zbar,
           snu, dsnudt, dsnudd, dsnuda, dsnudz);

    sp(i, j, k, vars.isneut) = snu;
    sp(i, j, k, vars.isneutdt) = dsnudt;
    sp(i, j, k, vars.isneutda) = dsnuda;
    sp(i, j, k, vars.isneutdz) = dsnudz;

  });

}
