#include <variables.H>
#include <network.H>

plot_t init_variables() {

  plot_t p;

  p.irho = p.next_index(1);
  p.itemp = p.next_index(1);
  p.ispec = p.next_index(NumSpec);
  p.ispec_old = p.next_index(NumSpec);
  p.irodot = p.next_index(NumSpec);
  p.irho_Hnuc = p.net_index(1);

  p.names.resize(p.n_plot_comps);

  p.names[p.irho] = "density";
  p.names[p.itemp] = "temperature";
  for (int n = 0; n < NumSpec; n++) {
    p.names[p.ispec + n] = "Xnew_" + spec_names_cxx[n];
    p.names[p.ispec_old + n] = "Xold_" + spec_names_cxx[n];
    p.names[p.irodot + n] = "wdot_" + spec_names_cxx[n];
  }
  p.names[p.irho_Huc] = "rho_Hnuc";

  return p;
}

