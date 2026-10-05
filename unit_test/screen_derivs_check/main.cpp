#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

using namespace amrex;

#include <screening_derivs.H>
#include <AMReX_buildInfo.H>

#include <network.H>
#include <eos.H>
#include <screen.H>

#include <cmath>
#include <unit_test.H>

using namespace unit_test_rp;

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

    // do the runtime parameter initializations and microphysics inits
    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "reading extern runtime parameters ..." << std::endl;
    }

    init_unit_test();

    eos_init(small_temp, small_dens);

    network_init();

    screening_init();

    test_screening_derivatives();
}


int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    main_main();

    amrex::Finalize();
    return 0;
}
