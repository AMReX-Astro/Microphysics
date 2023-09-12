#include <variables.H>
#include <network.H>

plot_t init_variables() {

  plot_t p;

  p.irho = p.next_index(1);
  p.itemp = p.next_index(1);
  p.ispec = p.next_index(NumSpec);
  p.isneut = p.next_index(1);
  p.isneutdt = p.next_index(1);
  p.isneutda = p.next_index(1);
  p.isneutdz = p.next_index(1);

  return p;
}


void get_varnames(const plot_t& p, amrex::Vector<std::string>& names) {

  names.resize(p.n_plot_comps);

  names[p.irho] = "density";
  names[p.itemp] = "temperature";
  for (int n = 0; n < NumSpec; n++) {
    names[p.ispec + n] = "X_" + spec_names_cxx[n];
  }

  names[p.isneut] = "sneut";
  names[p.isneutdt] = "dsneutdt";
  names[p.isneutda] = "dsneutda";
  names[p.isneutdz] = "dsneutdz";
}


