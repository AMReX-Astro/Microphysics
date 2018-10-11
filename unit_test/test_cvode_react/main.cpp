#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_BCRec.H>


using namespace amrex;

#include "test_react.H"
#include "test_react_F.H"
#include "AMReX_buildInfo.H"


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
    amrex::Real tmax;
    Vector<int> bc_lo(AMREX_SPACEDIM,0);
    Vector<int> bc_hi(AMREX_SPACEDIM,0);

    std::string prefix = "plt";

    IntVect tile_size(1024, 8, 8);

    // inputs parameters
    {
        // ParmParse is way of reading inputs from the inputs file
        ParmParse pp;

        // We need to get n_cell from the inputs file - this is the
        // number of cells on each side of a square (or cubic) domain.
        pp.get("n_cell", n_cell);

	// Get the integration time
	pp.get("tmax", tmax);

        // The domain is broken into boxes of size max_grid_size
        max_grid_size = 32;
        pp.query("max_grid_size", max_grid_size);

        pp.query("prefix", prefix);

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

    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "initializing unit test ..." << std::endl;
    }

    cudaError_t cuda_status = cudaSuccess;
    cuda_status = cudaDeviceSynchronize();
    assert(cuda_status == cudaSuccess);
    
    init_unit_test(probin_file_name.dataPtr(), &probin_file_length);

    cuda_status = cudaDeviceSynchronize();
    assert(cuda_status == cudaSuccess);    

    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "finished initializing unit test ..." << std::endl;
    }
    
    // Ncomp = number of components for each array
    int Ncomp = -1;
    init_variables();
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

    // time = starting time in the simulation
    Real time = 0.0;

    // How Boxes are distrubuted among MPI processes
    DistributionMapping dm(ba);

    // we allocate our main multifabs
    MultiFab state(ba, dm, Ncomp, Nghost);

    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "initializing state ..." << std::endl;
    }
    
    // Initialize the state
    for ( MFIter mfi(state); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.validbox();

        init_state(AMREX_ARLIM_ANYD(bx.loVect()), AMREX_ARLIM_ANYD(bx.hiVect()),
                   BL_TO_FORTRAN_ANYD(state[mfi]), &n_cell);

    }

    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "finished initializing state ..." << std::endl;
    }

    // What time is it now?  We'll use this to compute total react time.
    Real strt_time = ParallelDescriptor::second();

    // allocate a multifab for the number of RHS calls
    // so we can manually do the reductions (for GPU)
    iMultiFab integrator_n_rhs(ba, dm, 1, Nghost);

    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "reacting state with timestep " << tmax << " ..." << std::endl;
    }
    
    // Do the reactions
#ifdef _OPENMP
#pragma omp parallel
#endif
    for ( MFIter mfi(state, tile_size); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.tilebox();

        do_react(AMREX_ARLIM_3D(bx.loVect()), AMREX_ARLIM_3D(bx.hiVect()),
		 BL_TO_FORTRAN_ANYD(state[mfi]), Ncomp, tmax);

    }

    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "finished reacting state ..." << std::endl;
    }

    int n_rhs_min = integrator_n_rhs.min(0);
    int n_rhs_max = integrator_n_rhs.max(0);
    int n_rhs_sum = integrator_n_rhs.norm1();

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

    // Write a plotfile
    int n = 0;

    WriteSingleLevelPlotfile(prefix + name + integrator, state, varnames, geom, time, 0);


    // Call the timer again and compute the maximum difference between
    // the start time and stop time over all processors
    Real stop_time = ParallelDescriptor::second() - strt_time;
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::ReduceRealMax(stop_time, IOProc);

    // Tell the I/O Processor to write out the "run time"
    amrex::Print() << "Run time = " << stop_time << std::endl;

    // print statistics
    std::cout << "min number of rhs calls: " << n_rhs_min << std::endl;
    std::cout << "avg number of rhs calls: " << n_rhs_sum / (n_cell*n_cell*n_cell) << std::endl;
    std::cout << "max number of rhs calls: " << n_rhs_max << std::endl;

}
