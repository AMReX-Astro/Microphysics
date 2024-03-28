#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>


using namespace amrex;

#include <test_parameters.H>
#include <AMReX_buildInfo.H>

#include <extern_parameters.H>
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

    // do the runtime parameter initializations and microphysics inits
    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "reading extern runtime parameters ..." << std::endl;
    }

    ParmParse ppa("amr");

    auto params = init_unit_test();

    std::cout << "in C++" << std::endl;

    std::cout << "  eos_input_is_constant = " << eos_rp::eos_input_is_constant << " " << params.eos.eos_input_is_constant << std::endl;
    std::cout << "  test_string = " << unit_test_rp::test_string << " " << params.unit_test.test_string << std::endl;
    std::cout << "  dens_min = " << unit_test_rp::dens_min << " " << params.unit_test.dens_min << std::endl;
    std::cout << "  nonaka_file = " << integrator_rp::nonaka_file << " " << params.integrator.nonaka_file << std::endl;

}
