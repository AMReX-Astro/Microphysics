#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>


using namespace amrex;

#include "test_eos.H"
#include "test_eos_F.H"
#include "AMReX_buildInfo.H"

#include <network.H>
#include <eos.H>
#include <variables.H>

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    main_main();

    amrex::Finalize();
    return 0;
}

void main_main ()
{

    // AMREX_SPACEDIM: number of dimensions
    int n_cell, max_grid_size;
    Vector<int> bc_lo(AMREX_SPACEDIM,0);
    Vector<int> bc_hi(AMREX_SPACEDIM,0);

    // inputs parameters
    {
        // ParmParse is way of reading inputs from the inputs file
        ParmParse pp;

        // We need to get n_cell from the inputs file - this is the
        // number of cells on each side of a square (or cubic) domain.
        pp.get("n_cell", n_cell);

        // The domain is broken into boxes of size max_grid_size
        max_grid_size = 32;
        pp.query("max_grid_size", max_grid_size);

    }

    Vector<int> is_periodic(AMREX_SPACEDIM,0);
    for (int idim=0; idim < AMREX_SPACEDIM; ++idim) {
      is_periodic[idim] = 1;
    }

    // make BoxArray and Geometry
    BoxArray ba;
    Geometry geom;
    {
        IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
        IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));
        Box domain(dom_lo, dom_hi);

        // Initialize the boxarray "ba" from the single box "bx"
        ba.define(domain);

        // Break up boxarray "ba" into chunks no larger than
        // "max_grid_size" along a direction
        ba.maxSize(max_grid_size);

        // This defines the physical box, [0, 1] in each direction.
        RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                         {AMREX_D_DECL(1.0, 1.0, 1.0)});

        // This defines a Geometry object
        geom.define(domain, &real_box,
                    CoordSys::cartesian, is_periodic.data());
    }

    // Nghost = number of ghost cells for each array
    int Nghost = 0;

    // do the runtime parameter initializations and microphysics inits
    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "reading extern runtime parameters ..." << std::endl;
    }

    ParmParse ppa("amr");

    std::string probin_file = "probin";

    ppa.query("probin_file", probin_file);

    const int probin_file_length = probin_file.length();
    Vector<int> probin_file_name(probin_file_length);

    for (int i = 0; i < probin_file_length; i++)
      probin_file_name[i] = probin_file[i];

    init_unit_test(probin_file_name.dataPtr(), &probin_file_length);

    init_extern_parameters();

    eos_init();

    auto vars = init_variables();

    // time = starting time in the simulation
    Real time = 0.0;

    // How Boxes are distrubuted among MPI processes
    DistributionMapping dm(ba);

    // we allocate our main multifabs
    MultiFab state(ba, dm, vars.n_plot_comps, Nghost);

    // Initialize the state to zero; we will fill
    // it in below in do_eos.
    state.setVal(0.0);

    // What time is it now?  We'll use this to compute total run time.
    Real strt_time = ParallelDescriptor::second();

    const int ih1 = network_spec_index("hydrogen-1");
    const int ihe4 = network_spec_index("helium-4");

    Real dlogrho = 0.0e0_rt;
    Real dlogT   = 0.0e0_rt;
    Real dmetal  = 0.0e0_rt;

    if (n_cell > 1) {
        dlogrho = (log10(dens_max) - log10(dens_min))/(n_cell - 1);
        dlogT   = (log10(temp_max) - log10(temp_min))/(n_cell - 1);
        dmetal  = (metalicity_max  - 0.0)/(n_cell - 1);
    }

    // Initialize the state and compute the different thermodynamics
    // by inverting the EOS
    for ( MFIter mfi(state); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.validbox();

        Array4<Real> const sp = state.array(mfi);

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

          Real temp_zone = std::pow(10.0, log10(temp_min) + static_cast<Real>(j)*dlogT);
          eos_state.T = temp_zone;

          Real dens_zone = std::pow(10.0, log10(dens_min) + static_cast<Real>(i)*dlogrho);
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

    // Call the timer again and compute the maximum difference between
    // the start time and stop time over all processors
    Real stop_time = ParallelDescriptor::second() - strt_time;
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::ReduceRealMax(stop_time, IOProc);


    std::string name = "test_eos_C.";

    // Write a plotfile
    WriteSingleLevelPlotfile(name + eos_name, state, vars.names, geom, time, 0);

    // Tell the I/O Processor to write out the "run time"
    amrex::Print() << "Run time = " << stop_time << std::endl;

}
