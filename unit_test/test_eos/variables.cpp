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
  p.ispec = p.next_index(NumSpec);

  p.ierr_T_eos_rh = p.next_index(1);
  p.ierr_rho_eos_tp = p.next_index(1);
  p.ierr_T_eos_rp = p.next_index(1);
  p.ierr_T_eos_re = p.next_index(1);
  p.ierr_rho_eos_ps = p.next_index(1);
  p.ierr_T_eos_ps = p.next_index(1);
  p.ierr_rho_eos_ph = p.next_index(1);
  p.ierr_T_eos_ph = p.next_index(1);
  p.ierr_rho_eos_th = p.next_index(1);

  p.icv = p.next_index(1);
  p.icp = p.next_index(1);
  p.ine = p.next_index(1);
  p.inp = p.next_index(1);
  p.ieta = p.next_index(1);
  p.ipele = p.next_index(1);
  p.ippos = p.next_index(1);
  p.imu = p.next_index(1);
  p.imue = p.next_index(1);
  p.iye = p.next_index(1);
  p.idpdt = p.next_index(1);
  p.idpdr = p.next_index(1);
  p.idedt = p.next_index(1);
  p.idedr = p.next_index(1);
  p.idsdt = p.next_index(1);
  p.idsdr = p.next_index(1);
  p.idhdt = p.next_index(1);
  p.idhdr = p.next_index(1);
  p.idpdx = p.next_index(NumSpec);
  p.idedx = p.next_index(NumSpec);
  p.idhdx = p.next_index(NumSpec);
  p.igam1 = p.next_index(1);
  p.ics = p.next_index(1);
  p.iabar = p.next_index(1);
  p.izbar = p.next_index(1);
  p.idpda = p.next_index(1);
  p.idpdz = p.next_index(1);
  p.ideda = p.next_index(1);
  p.idedz = p.next_index(1);
  p.idpde = p.next_index(1);
  p.idpdre = p.next_index(1);

  return p;
}


void get_varnames(const plot_t p, amrex::Vector<std::string>& names) {

  names.resize(p.n_plot_comps);

  names[p.irho] = "density";
  names[p.itemp] = "temperature";
  names[p.ih] = "specific_enthalpy";
  names[p.ie] = "specific_energy";
  names[p.ip] = "pressure";
  names[p.is] = "specific_entropy";
  for (int n = 0; n < NumSpec; n++) {
    names[p.ispec + n] = "X_" + spec_names_cxx[n];
  }

  names[p.ierr_T_eos_rh] = "err_T_eos_rh";
  names[p.ierr_rho_eos_tp] = "err_rho_eos_tp";
  names[p.ierr_T_eos_rp] = "err_T_eos_rp";
  names[p.ierr_T_eos_re] = "err_T_eos_re";
  names[p.ierr_rho_eos_ps] = "err_rho_eos_ps";
  names[p.ierr_T_eos_ps] = "err_T_eos_ps";
  names[p.ierr_rho_eos_ph] = "err_rho_eos_ph";
  names[p.ierr_T_eos_ph] = "err_T_eos_ph";
  names[p.ierr_rho_eos_th] = "err_rho_eos_th";

  names[p.icv] = "c_v";
  names[p.icp] = "c_p";
  names[p.ine] = "n_e";
  names[p.inp] = "n_p";
  names[p.ieta] = "eta";
  names[p.ipele] = "p_ele";
  names[p.ippos] = "p_pos";
  names[p.imu] = "mu";
  names[p.imue] = "mu_e";
  names[p.iye] = "Y_e";
  names[p.idpdt] = "dp_dT";
  names[p.idpdr] = "dp_drho";
  names[p.idedt] = "de_dT";
  names[p.idedr] = "de_drho";
  names[p.idsdt] = "ds_dT";
  names[p.idsdr] = "ds_drho";
  names[p.idhdt] = "dh_dT";
  names[p.idhdr] = "dh_drho";
  for (int n = 0; n < NumSpec; n++) {
    names[p.idpdx + n] = "dp_dX_" + spec_names_cxx[n];
    names[p.idedx + n] = "de_dX_" + spec_names_cxx[n];
    names[p.idhdx + n] = "dh_dX_" + spec_names_cxx[n];
  }
  names[p.igam1] = "Gamma_1";
  names[p.ics] = "soundspeed";
  names[p.iabar] = "Abar";
  names[p.izbar] = "Zbar";
  names[p.idpda] = "dp_dA";
  names[p.idpdz] = "dp_dZ";
  names[p.ideda] = "de_dA";
  names[p.idedz] = "de_dZ";
  names[p.idpde] = "dp_de_rho";
  names[p.idpdre] = "dp_drho_e";

}


