#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>


using namespace amrex;

#include <test_parameters.H>
#include <test_parameters_F.H>
#include <AMReX_buildInfo.H>

#include <extern_parameters.H>
#include <unit_test_F.H>
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

    std::string probin_file = "probin";

    ppa.query("probin_file", probin_file);

    const int probin_file_length = probin_file.length();
    Vector<int> probin_file_name(probin_file_length);

    for (int i = 0; i < probin_file_length; i++)
      probin_file_name[i] = probin_file[i];

    // initialize the F90 parameters
    init_unit_test(probin_file_name.dataPtr(), &probin_file_length);

    do_f90_parameters();

    std::cout << "in C++" << std::endl;

    std::cout << "  eos_input_is_constant = " << eos_input_is_constant << std::endl;
    std::cout << "  test_string = " << test_string << std::endl;
    std::cout << "  dens_min = " << dens_min << std::endl;
    std::cout << "  nonaka_file = " << nonaka_file << std::endl;

}
