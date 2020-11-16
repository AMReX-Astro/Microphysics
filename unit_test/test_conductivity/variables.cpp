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


  return p;
}


void get_varnames(const plot_t p, amrex::Vector<std::string>& names) {

  names.resize(p.n_plot_comps);

  names[p.irho] = "density";
  names[p.itemp] = "temperature";
  names[p.ih] = "specific_enthalpy";
  names[p.ie] = "specific_energy";
  names[p.ip] = "pressure";
  names[p.iconductivity] = "conductivity";
  names[p.is] = "specific_entropy";
  for (int n = 0; n < NumSpec; n++) {
    names[p.ispec + n] = "X_" + spec_names_cxx[n];
  }

}
