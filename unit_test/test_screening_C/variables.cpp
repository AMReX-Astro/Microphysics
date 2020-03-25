#include <variables.H>
#include <network.H>

plot_t init_variables() {

  plot_t p;

  p.irho = p.next_index(1);
  p.itemp = p.next_index(1);
  p.ispec = p.next_index(NumSpec);

  p.iscn_he4_he4 = p.next_index(1);
  p.iscn_he4_be8 = p.next_index(1);
  p.iscn_c12_he4 = p.next_index(1);
  p.iscn_c12_c12 = p.next_index(1);
  p.iscn_c12_o16 = p.next_index(1);
  p.iscn_o16_o16 = p.next_index(1);
  p.iscn_o16_he4 = p.next_index(1);
  p.iscn_ne20_he4 = p.next_index(1);
  p.iscn_mg24_he4 = p.next_index(1);
  p.iscn_al27_p = p.next_index(1);
  p.iscn_si28_he4 = p.next_index(1);
  p.iscn_p31_p = p.next_index(1);
  p.iscn_s32_he4 = p.next_index(1);
  p.iscn_cl35_p = p.next_index(1);
  p.iscn_ar36_he4 = p.next_index(1);
  p.iscn_k39_p = p.next_index(1);
  p.iscn_ca40_he4 = p.next_index(1);
  p.iscn_sc43_p = p.next_index(1);
  p.iscn_ti44_he4 = p.next_index(1);
  p.iscn_v47_p = p.next_index(1);
  p.iscn_cr48_he4 = p.next_index(1);
  p.iscn_mn51_p = p.next_index(1);
  p.iscn_fe52_he4 = p.next_index(1);
  p.iscn_co55_p = p.next_index(1);
  p.iscn_fe54_p = p.next_index(1);
  p.iscn_fe54_he4 = p.next_index(1);
  p.iscn_fe56_p = p.next_index(1);
  p.iscn_d_p = p.next_index(1);
  p.iscn_p_p = p.next_index(1);
  p.iscn_he3_he3 = p.next_index(1);
  p.iscn_he3_he4 = p.next_index(1);
  p.iscn_c12_p = p.next_index(1);
  p.iscn_n14_p = p.next_index(1);
  p.iscn_o16_p = p.next_index(1);
  p.iscn_n14_he4 = p.next_index(1);

  p.iscn_he4_he4_dt = p.next_index(1);
  p.iscn_he4_be8_dt = p.next_index(1);
  p.iscn_c12_he4_dt = p.next_index(1);
  p.iscn_c12_c12_dt = p.next_index(1);
  p.iscn_c12_o16_dt = p.next_index(1);
  p.iscn_o16_o16_dt = p.next_index(1);
  p.iscn_o16_he4_dt = p.next_index(1);
  p.iscn_ne20_he4_dt = p.next_index(1);
  p.iscn_mg24_he4_dt = p.next_index(1);
  p.iscn_al27_p_dt = p.next_index(1);
  p.iscn_si28_he4_dt = p.next_index(1);
  p.iscn_p31_p_dt = p.next_index(1);
  p.iscn_s32_he4_dt = p.next_index(1);
  p.iscn_cl35_p_dt = p.next_index(1);
  p.iscn_ar36_he4_dt = p.next_index(1);
  p.iscn_k39_p_dt = p.next_index(1);
  p.iscn_ca40_he4_dt = p.next_index(1);
  p.iscn_sc43_p_dt = p.next_index(1);
  p.iscn_ti44_he4_dt = p.next_index(1);
  p.iscn_v47_p_dt = p.next_index(1);
  p.iscn_cr48_he4_dt = p.next_index(1);
  p.iscn_mn51_p_dt = p.next_index(1);
  p.iscn_fe52_he4_dt = p.next_index(1);
  p.iscn_co55_p_dt = p.next_index(1);
  p.iscn_fe54_p_dt = p.next_index(1);
  p.iscn_fe54_he4_dt = p.next_index(1);
  p.iscn_fe56_p_dt = p.next_index(1);
  p.iscn_d_p_dt = p.next_index(1);
  p.iscn_p_p_dt = p.next_index(1);
  p.iscn_he3_he3_dt = p.next_index(1);
  p.iscn_he3_he4_dt = p.next_index(1);
  p.iscn_c12_p_dt = p.next_index(1);
  p.iscn_n14_p_dt = p.next_index(1);
  p.iscn_o16_p_dt = p.next_index(1);
  p.iscn_n14_he4_dt = p.next_index(1);

  p.names.resize(p.n_plot_comps);

  p.names[p.irho] = "density";
  p.names[p.itemp] = "temperature";
  for (int n = 0; n < NumSpec; n++) {
    p.names[p.ispec + n] = "X_" + spec_names_cxx[n];
  }

  p.names(p.iscn_he4_he4) = "scn_he4_he4";
  p.names(p.iscn_he4_be8) = "scn_he4_be8";
  p.names(p.iscn_c12_he4) = "scn_c12_he4";
  p.names(p.iscn_c12_c12) = "scn_c12_c12";
  p.names(p.iscn_c12_o16) = "scn_c12_o16";
  p.names(p.iscn_o16_o16) = "scn_o16_o16";
  p.names(p.iscn_o16_he4) = "scn_o16_he4";
  p.names(p.iscn_ne20_he4) = "scn_ne20_he4";
  p.names(p.iscn_mg24_he4) = "scn_mg24_he4";
  p.names(p.iscn_al27_p) = "scn_al27_p";
  p.names(p.iscn_si28_he4) = "scn_si28_he4";
  p.names(p.iscn_p31_p) = "scn_p31_p";
  p.names(p.iscn_s32_he4) = "scn_s32_he4";
  p.names(p.iscn_cl35_p) = "scn_cl35_p";
  p.names(p.iscn_ar36_he4) = "scn_ar36_he4";
  p.names(p.iscn_k39_p) = "scn_k39_p";
  p.names(p.iscn_ca40_he4) = "scn_ca40_he4";
  p.names(p.iscn_sc43_p) = "scn_sc43_p";
  p.names(p.iscn_ti44_he4) = "scn_ti44_he4";
  p.names(p.iscn_v47_p) = "scn_v47_p";
  p.names(p.iscn_cr48_he4) = "scn_cr48_he4";
  p.names(p.iscn_mn51_p) = "scn_mn51_p";
  p.names(p.iscn_fe52_he4) = "scn_fe52_he4";
  p.names(p.iscn_co55_p) = "scn_co55_p";
  p.names(p.iscn_fe54_p) = "scn_fe54_p";
  p.names(p.iscn_fe54_he4) = "scn_fe54_he4";
  p.names(p.iscn_fe56_p) = "scn_fe56_p";
  p.names(p.iscn_d_p) = "scn_d_p";
  p.names(p.iscn_p_p) = "scn_p_p";
  p.names(p.iscn_he3_he3) = "scn_he3_he3";
  p.names(p.iscn_he3_he4) = "scn_he3_he4";
  p.names(p.iscn_c12_p) = "scn_c12_p";
  p.names(p.iscn_n14_p) = "scn_n14_p";
  p.names(p.iscn_o16_p) = "scn_o16_p";
  p.names(p.iscn_n14_he4) = "scn_n14_he4";

  p.names(p.iscn_he4_he4_dt) = "scn_he4_he4_dt";
  p.names(p.iscn_he4_be8_dt) = "scn_he4_be8_dt";
  p.names(p.iscn_c12_he4_dt) = "scn_c12_he4_dt";
  p.names(p.iscn_c12_c12_dt) = "scn_c12_c12_dt";
  p.names(p.iscn_c12_o16_dt) = "scn_c12_o16_dt";
  p.names(p.iscn_o16_o16_dt) = "scn_o16_o16_dt";
  p.names(p.iscn_o16_he4_dt) = "scn_o16_he4_dt";
  p.names(p.iscn_ne20_he4_dt) = "scn_ne20_he4_dt";
  p.names(p.iscn_mg24_he4_dt) = "scn_mg24_he4_dt";
  p.names(p.iscn_al27_p_dt) = "scn_al27_p_dt";
  p.names(p.iscn_si28_he4_dt) = "scn_si28_he4_dt";
  p.names(p.iscn_p31_p_dt) = "scn_p31_p_dt";
  p.names(p.iscn_s32_he4_dt) = "scn_s32_he4_dt";
  p.names(p.iscn_cl35_p_dt) = "scn_cl35_p_dt";
  p.names(p.iscn_ar36_he4_dt) = "scn_ar36_he4_dt";
  p.names(p.iscn_k39_p_dt) = "scn_k39_p_dt";
  p.names(p.iscn_ca40_he4_dt) = "scn_ca40_he4_dt";
  p.names(p.iscn_sc43_p_dt) = "scn_sc43_p_dt";
  p.names(p.iscn_ti44_he4_dt) = "scn_ti44_he4_dt";
  p.names(p.iscn_v47_p_dt) = "scn_v47_p_dt";
  p.names(p.iscn_cr48_he4_dt) = "scn_cr48_he4_dt";
  p.names(p.iscn_mn51_p_dt) = "scn_mn51_p_dt";
  p.names(p.iscn_fe52_he4_dt) = "scn_fe52_he4_dt";
  p.names(p.iscn_co55_p_dt) = "scn_co55_p_dt";
  p.names(p.iscn_fe54_p_dt) = "scn_fe54_p_dt";
  p.names(p.iscn_fe54_he4_dt) = "scn_fe54_he4_dt";
  p.names(p.iscn_fe56_p_dt) = "scn_fe56_p_dt";
  p.names(p.iscn_d_p_dt) = "scn_d_p_dt";
  p.names(p.iscn_p_p_dt) = "scn_p_p_dt";
  p.names(p.iscn_he3_he3_dt) = "scn_he3_he3_dt";
  p.names(p.iscn_he3_he4_dt) = "scn_he3_he4_dt";
  p.names(p.iscn_c12_p_dt) = "scn_c12_p_dt";
  p.names(p.iscn_n14_p_dt) = "scn_n14_p_dt";
  p.names(p.iscn_o16_p_dt) = "scn_o16_p_dt";
  p.names(p.iscn_n14_he4_dt) = "scn_n14_he4_dt";

  return p;
}

