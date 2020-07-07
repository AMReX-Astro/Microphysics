#include <variables.H>
#include <network.H>

plot_t init_variables() {

  plot_t p;

  p.irho = p.next_index(1);
  p.itemp = p.next_index(1);
  p.ih = p.next_index(1);
  p.ie = p.next_index(1);
  p.ip = p.next_index(1);
  p.is = p.next_index(1);
  p.iconductivity = p.next_index(1);
  p.ispec = p.next_index(NumSpec);

  p.names.resize(p.n_plot_comps);

  p.names[p.irho] = "density";
  p.names[p.itemp] = "temperature";
  p.names[p.ih] = "specific_enthalpy";
  p.names[p.ie] = "specific_energy";
  p.names[p.ip] = "pressure";
  p.names[p.iconductivity] = "conductivity";
  p.names[p.is] = "specific_entropy";
  for (int n = 0; n < NumSpec; n++) {
    p.names[p.ispec + n] = "X_" + spec_names_cxx[n];
  }

  return p;
}

