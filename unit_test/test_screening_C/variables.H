#include <vector>
#include <string>

#ifndef _variables_H_
#define _variables_H_

#include <AMReX_Vector.H>

class plot_t {

public:


  int irho = -1;
  int itemp = -1;
  int ispec = -1;

  int iscn_he4_he4 = -1;
  int iscn_he4_be8 = -1;
  int iscn_c12_he4 = -1;
  int iscn_c12_c12 = -1;
  int iscn_c12_o16 = -1;
  int iscn_o16_o16 = -1;
  int iscn_o16_he4 = -1;
  int iscn_ne20_he4 = -1;
  int iscn_mg24_he4 = -1;
  int iscn_al27_p = -1;
  int iscn_si28_he4 = -1;
  int iscn_p31_p = -1;
  int iscn_s32_he4 = -1;
  int iscn_cl35_p = -1;
  int iscn_ar36_he4 = -1;
  int iscn_k39_p = -1;
  int iscn_ca40_he4 = -1;
  int iscn_sc43_p = -1;
  int iscn_ti44_he4 = -1;
  int iscn_v47_p = -1;
  int iscn_cr48_he4 = -1;
  int iscn_mn51_p = -1;
  int iscn_fe52_he4 = -1;
  int iscn_co55_p = -1;
  int iscn_fe54_p = -1;
  int iscn_fe54_he4 = -1;
  int iscn_fe56_p = -1;
  int iscn_d_p = -1;
  int iscn_p_p = -1;
  int iscn_he3_he3 = -1;
  int iscn_he3_he4 = -1;
  int iscn_c12_p = -1;
  int iscn_n14_p = -1;
  int iscn_o16_p = -1;
  int iscn_n14_he4 = -1;

  int iscn_he4_he4_dt = -1;
  int iscn_he4_be8_dt = -1;
  int iscn_c12_he4_dt = -1;
  int iscn_c12_c12_dt = -1;
  int iscn_c12_o16_dt = -1;
  int iscn_o16_o16_dt = -1;
  int iscn_o16_he4_dt = -1;
  int iscn_ne20_he4_dt = -1;
  int iscn_mg24_he4_dt = -1;
  int iscn_al27_p_dt = -1;
  int iscn_si28_he4_dt = -1;
  int iscn_p31_p_dt = -1;
  int iscn_s32_he4_dt = -1;
  int iscn_cl35_p_dt = -1;
  int iscn_ar36_he4_dt = -1;
  int iscn_k39_p_dt = -1;
  int iscn_ca40_he4_dt = -1;
  int iscn_sc43_p_dt = -1;
  int iscn_ti44_he4_dt = -1;
  int iscn_v47_p_dt = -1;
  int iscn_cr48_he4_dt = -1;
  int iscn_mn51_p_dt = -1;
  int iscn_fe52_he4_dt = -1;
  int iscn_co55_p_dt = -1;
  int iscn_fe54_p_dt = -1;
  int iscn_fe54_he4_dt = -1;
  int iscn_fe56_p_dt = -1;
  int iscn_d_p_dt = -1;
  int iscn_p_p_dt = -1;
  int iscn_he3_he3_dt = -1;
  int iscn_he3_he4_dt = -1;
  int iscn_c12_p_dt = -1;
  int iscn_n14_p_dt = -1;
  int iscn_o16_p_dt = -1;
  int iscn_n14_he4_dt = -1;

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
