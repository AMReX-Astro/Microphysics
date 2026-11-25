#include <iostream>

#include <extern_parameters.H>
#include <eos.H>
#include <network.H>
#include <plot_eos.H>

int main(int argc, char *argv[]) {

  amrex::Initialize(argc, argv);

  // initialize runtime parameters

  init_extern_parameters();

  // initialize EOS (must be done after init_extern_parameters)

  eos_init();

  // initialize network (this may set nuclei properties needed by the EOS)

  network_init();

  output_plot_info();

  amrex::Finalize();
}
