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

#include <eos.H>
#include <variables.H>
#include <actual_eos.H>

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

    eos_cxx_init();

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

    Real dlogrho = (log10(dens_max) - log10(dens_min))/(n_cell - 1);
    Real dlogT = (log10(temp_max) - log10(temp_min))/(n_cell - 1);
    Real dmetal = (metalicity_max  - 0.0)/(n_cell - 1);

    // Initialize the state and compute the different thermodynamics
    // by inverting the EOS
    for ( MFIter mfi(state); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.validbox();

        Array4<Real> const state_arr = state.array(mfi);

        amrex::Print() << "here!!!!" << std::endl;

        AMREX_PARALLEL_FOR_3D(bx, i, j, k,
        {

          // set the composition -- approximately solar
          Real metalicity = 0.0 + static_cast<Real> (k) * dmetal;

          eos_t eos_state;

          for (int n = 0; n < NumSpec; n++) {
            eos_state.xn[n] = metalicity/(NumSpec - 2);
          }
          eos_state.xn[0] = 0.75 - 0.5*metalicity;
          eos_state.xn[1] = 0.25 - 0.5*metalicity;

          eos_state.T = std::pow(10.0, log10(temp_min) + static_cast<double>(j)*dlogT);
          eos_state.rho = std::pow(10.0, log10(dens_min) + static_cast<double>(i)*dlogrho);

          eos_cxx(eos_input_rt, eos_state);

          amrex::Print() << eos_state.p << std::endl;
        });



    }

    // Call the timer again and compute the maximum difference between
    // the start time and stop time over all processors
    Real stop_time = ParallelDescriptor::second() - strt_time;
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::ReduceRealMax(stop_time, IOProc);


    std::string name = "test_eos.";

    // Write a plotfile
    WriteSingleLevelPlotfile(name + eos_name, state, vars.names, geom, time, 0);

    // Tell the I/O Processor to write out the "run time"
    amrex::Print() << "Run time = " << stop_time << std::endl;

}
