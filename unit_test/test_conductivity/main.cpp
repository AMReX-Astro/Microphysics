#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>


using namespace amrex;

#include <test_cond.H>
#include <AMReX_buildInfo.H>

#include <network.H>
#include <eos.H>
#include <conductivity.H>
#include <variables.H>

#include <cmath>

#include <unit_test.H>

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

    init_unit_test();

    eos_init(unit_test_rp::small_temp, unit_test_rp::small_dens);
    network_init();
    conductivity_init();

    // for C++
    plot_t vars;

    // C++ test
    vars = init_variables();

    amrex::Vector<std::string> names;
    get_varnames(vars, names);

    // time = starting time in the simulation
    Real time = 0.0;

    // How Boxes are distributed among MPI processes
    DistributionMapping dm(ba);

    // we allocate our main multifabs
    MultiFab state(ba, dm, vars.n_plot_comps, Nghost);

    // Initialize the state to zero; we will fill
    // it in below in do_eos.
    state.setVal(0.0);

    // What time is it now?  We'll use this to compute total run time.
    Real strt_time = ParallelDescriptor::second();


    Real dlogrho = 0.0e0_rt;
    Real dlogT   = 0.0e0_rt;
    Real dmetal  = 0.0e0_rt;

    if (n_cell > 1) {
        dlogrho = (std::log10(unit_test_rp::dens_max) - std::log10(unit_test_rp::dens_min))/(n_cell - 1);
        dlogT   = (std::log10(unit_test_rp::temp_max) - std::log10(unit_test_rp::temp_min))/(n_cell - 1);
        dmetal  = (unit_test_rp::metalicity_max  - 0.0)/(n_cell - 1);
    }

    // Initialize the state and compute the different thermodynamics
    // by inverting the EOS
    for ( MFIter mfi(state); mfi.isValid(); ++mfi )
    {
      const Box& bx = mfi.validbox();

      Array4<Real> const sp = state.array(mfi);

      cond_test_C(bx, dlogrho, dlogT, dmetal, vars, sp);

    }

    // Call the timer again and compute the maximum difference between
    // the start time and stop time over all processors
    Real stop_time = ParallelDescriptor::second() - strt_time;
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::ReduceRealMax(stop_time, IOProc);

    std::string name = "test_conductivity.";

    // Write a plotfile
    WriteSingleLevelPlotfile(name + std::string(cond_name),
                             state, names, geom, time, 0);

    write_job_info(name + std::string(cond_name));

    // Tell the I/O Processor to write out the "run time"
    amrex::Print() << "Run time = " << stop_time << std::endl;

}
