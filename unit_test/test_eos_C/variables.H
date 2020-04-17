#include <vector>
#include <string>

#ifndef _variables_H_
#define _variables_H_

#include <AMReX_Vector.H>

class plot_t {

public:


  int irho = -1;
  int itemp = -1;
  int ih = -1;
  int ie = -1;
  int ip = -1;
  int is = -1;
  int ispec = -1;

  int ierr_T_eos_rh = -1;
  int ierr_rho_eos_tp = -1;
  int ierr_T_eos_rp = -1;
  int ierr_T_eos_re = -1;
  int ierr_rho_eos_ps = -1;
  int ierr_T_eos_ps = -1;
  int ierr_rho_eos_ph = -1;
  int ierr_T_eos_ph = -1;
  int ierr_rho_eos_th = -1;

  int icv = -1;
  int icp = -1;
  int ine = -1;
  int inp = -1;
  int ieta = -1;
  int ipele = -1;
  int ippos = -1;
  int imu = -1;
  int imue = -1;
  int iye = -1;
  int idpdt = -1;
  int idpdr = -1;
  int idedt = -1;
  int idedr = -1;
  int idsdt = -1;
  int idsdr = -1;
  int idhdt = -1;
  int idhdr = -1;
  int idpdx = -1;
  int idedx = -1;
  int idhdx = -1;
  int igam1 = -1;
  int ics = -1;
  int iabar = -1;
  int izbar = -1;
  int idpda = -1;
  int idpdz = -1;
  int ideda = -1;
  int idedz = -1;
  int idpde = -1;
  int idpdre = -1;

  int n_plot_comps = 0;

  amrex::Vector<std::string> names;

  int next_index(const int num) {
    int next = n_plot_comps;
    n_plot_comps += num;
    return next;
  }

};

plot_t init_variables();


#endif
