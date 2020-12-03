#include <variables.H>
#include <network.H>

plot_t init_variables() {

  plot_t p;

  p.irho = p.next_index(1);
  p.itemp = p.next_index(1);

  p.ispec = p.next_index(NumSpec);
  p.ispec_old = p.next_index(NumSpec);

  p.itemp_dot = p.next_index(1);
  p.ienuc_dot = p.next_index(1);

  p.ijac = p.next_index(neqs * neqs);

  return p;

}


void get_varnames(const plot_t p, amrex::Vector<std::string>& names) {

  names.resize(p.n_plot_comps);

  names[p.irho] = "density";
  names[p.itemp] = "temperature";
  for (int n = 0; n < NumSpec; n++) {
    names[p.ispec + n] = "Ydot_" + spec_names_cxx[n];
    names[p.ispec_old + n] = "Xold_" + spec_names_cxx[n];
  }

  names[p.itemp_dot] = "Tdot";
  names[p.ienuc_dot] = "Edot";

  int n = 0;
  for (int j = 0; j < neqs; j++) {
      for (int i = 0; i < neqs; i++) {

          names[p.ijac + n] = "J_";

          if (i < NumSpec) {
              names[p.ijac + n] += spec_names_cxx[i] + "_";
          } else if (i == NumSpec) {
              names[p.ijac + n] += "T_";
          } else {
              names[p.ijac + n] += "E_";
          }

          if (j < NumSpec) {
              names[p.ijac + n] += spec_names_cxx[j];
          } else if (j == NumSpec) {
              names[p.ijac + n] += "T";
          } else {
              names[p.ijac + n] += "E";
          }

          n++;
      }
  }

}

