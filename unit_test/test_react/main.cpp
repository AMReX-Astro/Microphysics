#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_BCRec.H>


using namespace amrex;

#include <test_react.H>
#include <test_react_F.H>
#include <extern_parameters.H>
#include <eos.H>
#include <network.H>
#ifdef CXX_REACTIONS
#include <react_zones.H>
#endif
#include <AMReX_buildInfo.H>
#include <variables.H>
#include <unit_test.H>

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        main_main();
    }
    amrex::Finalize();
    return 0;
}

void main_main ()
{

    // AMREX_SPACEDIM: number of dimensions
    int n_cell, max_grid_size, print_every_nrhs;
    Vector<int> bc_lo(AMREX_SPACEDIM,0);
    Vector<int> bc_hi(AMREX_SPACEDIM,0);

    std::string prefix = "plt";

    IntVect tile_size(1024, 8, 8);

    int do_cxx = 0;

    // inputs parameters
    {
        // ParmParse is way of reading inputs from the inputs file
        ParmParse pp;

        // We need to get n_cell from the inputs file - this is the
        // number of cells on each side of a square (or cubic) domain.
        pp.get("n_cell", n_cell);

	print_every_nrhs = 0;
	pp.query("print_every_nrhs", print_every_nrhs);

        // The domain is broken into boxes of size max_grid_size
        max_grid_size = 32;
        pp.query("max_grid_size", max_grid_size);

        pp.query("prefix", prefix);

#ifdef CXX_REACTIONS
        pp.query("do_cxx", do_cxx);
#endif

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

    // Copy extern parameters from Fortran to C++
    init_extern_parameters();

    // C++ EOS initialization (must be done after Fortran eos_init and init_extern_parameters)
    eos_init();

#ifdef CXX_REACTIONS
    // C++ Network, RHS, screening, rates initialization
    network_init();
#endif

    // Ncomp = number of components for each array
    int Ncomp = -1;
    init_variables_F();
    get_ncomp(&Ncomp);

    int name_len = -1;
    get_name_len(&name_len);

    // get the variable names
    Vector<std::string> varnames;

    for (int i=0; i<Ncomp; i++) {
      char* cstring[name_len+1];
      get_var_name(cstring, &i);
      std::string name(*cstring);
      varnames.push_back(name);
    }

    plot_t vars;
    if (do_cxx == 1) {
      vars = init_variables();
    }

    // time = starting time in the simulation
    Real time = 0.0;

    // How Boxes are distrubuted among MPI processes
    DistributionMapping dm(ba);

    // we allocate our main multifabs
    MultiFab state(ba, dm, Ncomp, Nghost);

    // Initialize the state
    for ( MFIter mfi(state); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.validbox();

        init_state(AMREX_ARLIM_ANYD(bx.loVect()), AMREX_ARLIM_ANYD(bx.hiVect()),
                   BL_TO_FORTRAN_ANYD(state[mfi]), &n_cell);

    }

    // allocate a multifab for the number of RHS calls
    // so we can manually do the reductions (for GPU)
    iMultiFab integrator_n_rhs(ba, dm, 1, Nghost);

    // What time is it now?  We'll use this to compute total react time.
    Real strt_time = ParallelDescriptor::second();

    int num_failed = 0;

    // Do the reactions
#ifdef _OPENMP
#pragma omp parallel
#endif
    for ( MFIter mfi(state, tile_size); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.tilebox();

#ifdef CXX_REACTIONS
        if (do_cxx) {

            auto s = state.array(mfi);
            auto n_rhs = integrator_n_rhs.array(mfi);

            int* num_failed_d = AMREX_MFITER_REDUCE_SUM(&num_failed);

            AMREX_PARALLEL_FOR_3D(bx, i, j, k,
            {
                bool success = do_react(i, j, k, s, n_rhs, vars);

                if (!success) {
                    Gpu::Atomic::Add(num_failed_d, 1);
                }
            });

        }
        else {
#endif

#pragma gpu
            do_react(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                     BL_TO_FORTRAN_ANYD(state[mfi]),
                     BL_TO_FORTRAN_ANYD(integrator_n_rhs[mfi]));

#ifdef CXX_REACTIONS
        }
#endif

	if (print_every_nrhs != 0)
	  print_nrhs(AMREX_ARLIM_ANYD(bx.loVect()), AMREX_ARLIM_ANYD(bx.hiVect()),
		     BL_TO_FORTRAN_ANYD(integrator_n_rhs[mfi]));
    }

    if (num_failed > 0) {
        amrex::Abort("Integration failed");
    }

    // Call the timer again and compute the maximum difference between
    // the start time and stop time over all processors
    Real stop_time = ParallelDescriptor::second() - strt_time;
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::ReduceRealMax(stop_time, IOProc);

    int n_rhs_min = integrator_n_rhs.min(0);
    int n_rhs_max = integrator_n_rhs.max(0);
    long n_rhs_sum = integrator_n_rhs.sum(0, 0, true);

    // get the name of the integrator from the build info functions
    // written at compile time.  We will append the name of the
    // integrator to the output file name
    int nmodules = buildInfoGetNumModules();

    int int_idx = -1;
    std::string INT_KEY = "INTEGRATOR";
    for (int i=1; i<=nmodules; i++) {
      if (INT_KEY == buildInfoGetModuleName(i)) {
        int_idx = i;
      }
    }

    std::string name = "test_react.";
    std::string integrator = buildInfoGetModuleVal(int_idx);

#ifdef CXX_REACTIONS
    std::string language = do_cxx == 1 ? ".cxx" : "";
#else
    std::string language = "";
#endif

    // Write a plotfile
    WriteSingleLevelPlotfile(prefix + name + integrator + language, state, varnames, geom, time, 0);

    write_job_info(prefix + name + integrator + language);

    // Tell the I/O Processor to write out the "run time"
    amrex::Print() << "Run time = " << stop_time << std::endl;

    // print statistics
    std::cout << "min number of rhs calls: " << n_rhs_min << std::endl;
    std::cout << "avg number of rhs calls: " << n_rhs_sum / (n_cell*n_cell*n_cell) << std::endl;
    std::cout << "max number of rhs calls: " << n_rhs_max << std::endl;

}
