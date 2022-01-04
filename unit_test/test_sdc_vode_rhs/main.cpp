#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_BCRec.H>


using namespace amrex;

#include <test_react.H>
#include <extern_parameters.H>
#include <eos.H>
#include <network.H>
#include <vode_rhs_test.H>
#include <AMReX_buildInfo.H>
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


    // do the runtime parameter initializations and microphysics inits
    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "reading extern runtime parameters ..." << std::endl;
    }

    ParmParse ppa("amr");

    init_unit_test();

    // C++ EOS initialization (must be done after Fortran eos_init and init_extern_parameters)
    eos_init(small_temp, small_dens);

    // C++ Network, RHS, screening, rates initialization
    network_init();

    do_vode_rhs();


}
