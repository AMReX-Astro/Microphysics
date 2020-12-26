#include <variables.H>
#include <network.H>

plot_t init_variables() {

  plot_t p;

  p.irho = p.next_index(1);
  p.itemp = p.next_index(1);

  p.ispec = p.next_index(NumSpec);
  p.ispec_old = p.next_index(NumSpec);

#if NAUX_NET > 0
  p.iaux = p.next_index(NumAux);
  p.iaux_old = p.next_index(NumAux);
#endif

  p.irodot = p.next_index(NumSpec);
  p.irho_hnuc = p.next_index(1);

  return p;

}


void get_varnames(const plot_t p, amrex::Vector<std::string>& names) {

  names.resize(p.n_plot_comps);

  names[p.irho] = "density";
  names[p.itemp] = "temperature";
  for (int n = 0; n < NumSpec; n++) {
    names[p.ispec + n] = "Xnew_" + spec_names_cxx[n];
    names[p.ispec_old + n] = "Xold_" + spec_names_cxx[n];
    names[p.irodot + n] = "wdot_" + spec_names_cxx[n];
  }

#if NAUX_NET > 0
  for (int n = 0; n < NumAux; n++) {
      names[p.iaux + n] = "Xnew_" + aux_names_cxx[n];
      names[p.iaux_old + n] = "Xold_" + aux_names_cxx[n];
  }
#endif

  names[p.irho_hnuc] = "rho_Hnuc";
}

