#include <extern_parameters.H>
#include <unit_test.H>
#include <test_autodiff_arrays.H>

using namespace unit_test_rp;

int main(int argc, char *argv[]) {

  amrex::Initialize(argc, argv);

  init_unit_test();

  test_autodiff_arrays();

  amrex::Finalize();
}
