#include <AMReX_Array.H>
#include <string>
#include <table_rates.H>
#include <AMReX_Print.H>

namespace rate_tables
{

    AMREX_GPU_MANAGED table_t j_Na21_Ne21_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Na21_Ne21_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Na21_Ne21_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Na21_Ne21_temp;

    AMREX_GPU_MANAGED table_t j_Ne21_Na21_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ne21_Na21_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ne21_Na21_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ne21_Na21_temp;

    AMREX_GPU_MANAGED table_t j_Na22_Ne22_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Na22_Ne22_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Na22_Ne22_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Na22_Ne22_temp;

    AMREX_GPU_MANAGED table_t j_Ne22_Na22_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ne22_Na22_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ne22_Na22_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ne22_Na22_temp;

    AMREX_GPU_MANAGED table_t j_Mg23_Na23_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Mg23_Na23_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Mg23_Na23_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Mg23_Na23_temp;

    AMREX_GPU_MANAGED table_t j_Na23_Mg23_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Na23_Mg23_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Na23_Mg23_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Na23_Mg23_temp;

    AMREX_GPU_MANAGED table_t j_Al25_Mg25_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Al25_Mg25_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Al25_Mg25_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Al25_Mg25_temp;

    AMREX_GPU_MANAGED table_t j_Mg25_Al25_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Mg25_Al25_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Mg25_Al25_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Mg25_Al25_temp;

    AMREX_GPU_MANAGED table_t j_Al26_Mg26_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Al26_Mg26_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Al26_Mg26_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Al26_Mg26_temp;

    AMREX_GPU_MANAGED table_t j_Mg26_Al26_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Mg26_Al26_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Mg26_Al26_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Mg26_Al26_temp;

    AMREX_GPU_MANAGED table_t j_P29_Si29_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_P29_Si29_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_P29_Si29_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_P29_Si29_temp;

    AMREX_GPU_MANAGED table_t j_Si29_P29_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Si29_P29_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Si29_P29_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Si29_P29_temp;

    AMREX_GPU_MANAGED table_t j_P30_Si30_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_P30_Si30_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_P30_Si30_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_P30_Si30_temp;

    AMREX_GPU_MANAGED table_t j_Si30_P30_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Si30_P30_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Si30_P30_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Si30_P30_temp;

    AMREX_GPU_MANAGED table_t j_P31_Si31_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_P31_Si31_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_P31_Si31_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_P31_Si31_temp;

    AMREX_GPU_MANAGED table_t j_Si31_P31_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Si31_P31_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Si31_P31_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Si31_P31_temp;

    AMREX_GPU_MANAGED table_t j_P32_S32_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_P32_S32_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_P32_S32_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_P32_S32_temp;

    AMREX_GPU_MANAGED table_t j_P32_Si32_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_P32_Si32_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_P32_Si32_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_P32_Si32_temp;

    AMREX_GPU_MANAGED table_t j_S32_P32_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_S32_P32_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_S32_P32_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_S32_P32_temp;

    AMREX_GPU_MANAGED table_t j_Si32_P32_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Si32_P32_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Si32_P32_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Si32_P32_temp;

    AMREX_GPU_MANAGED table_t j_Cl33_S33_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Cl33_S33_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Cl33_S33_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Cl33_S33_temp;

    AMREX_GPU_MANAGED table_t j_P33_S33_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_P33_S33_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_P33_S33_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_P33_S33_temp;

    AMREX_GPU_MANAGED table_t j_S33_Cl33_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_S33_Cl33_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_S33_Cl33_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_S33_Cl33_temp;

    AMREX_GPU_MANAGED table_t j_S33_P33_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_S33_P33_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_S33_P33_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_S33_P33_temp;

    AMREX_GPU_MANAGED table_t j_Cl34_S34_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Cl34_S34_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Cl34_S34_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Cl34_S34_temp;

    AMREX_GPU_MANAGED table_t j_S34_Cl34_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_S34_Cl34_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_S34_Cl34_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_S34_Cl34_temp;

    AMREX_GPU_MANAGED table_t j_Cl35_S35_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Cl35_S35_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Cl35_S35_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Cl35_S35_temp;

    AMREX_GPU_MANAGED table_t j_S35_Cl35_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_S35_Cl35_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_S35_Cl35_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_S35_Cl35_temp;

    AMREX_GPU_MANAGED table_t j_Ar36_Cl36_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ar36_Cl36_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ar36_Cl36_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ar36_Cl36_temp;

    AMREX_GPU_MANAGED table_t j_Cl36_Ar36_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Cl36_Ar36_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Cl36_Ar36_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Cl36_Ar36_temp;

    AMREX_GPU_MANAGED table_t j_Cl36_S36_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Cl36_S36_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Cl36_S36_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Cl36_S36_temp;

    AMREX_GPU_MANAGED table_t j_S36_Cl36_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_S36_Cl36_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_S36_Cl36_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_S36_Cl36_temp;

    AMREX_GPU_MANAGED table_t j_Ar37_Cl37_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ar37_Cl37_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ar37_Cl37_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ar37_Cl37_temp;

    AMREX_GPU_MANAGED table_t j_Ar37_K37_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ar37_K37_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ar37_K37_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ar37_K37_temp;

    AMREX_GPU_MANAGED table_t j_Cl37_Ar37_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Cl37_Ar37_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Cl37_Ar37_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Cl37_Ar37_temp;

    AMREX_GPU_MANAGED table_t j_K37_Ar37_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_K37_Ar37_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_K37_Ar37_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_K37_Ar37_temp;

    AMREX_GPU_MANAGED table_t j_Ar38_K38_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ar38_K38_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ar38_K38_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ar38_K38_temp;

    AMREX_GPU_MANAGED table_t j_K38_Ar38_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_K38_Ar38_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_K38_Ar38_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_K38_Ar38_temp;

    AMREX_GPU_MANAGED table_t j_Ar39_K39_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ar39_K39_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ar39_K39_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ar39_K39_temp;

    AMREX_GPU_MANAGED table_t j_K39_Ar39_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_K39_Ar39_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_K39_Ar39_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_K39_Ar39_temp;

    AMREX_GPU_MANAGED table_t j_Ar40_K40_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ar40_K40_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ar40_K40_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ar40_K40_temp;

    AMREX_GPU_MANAGED table_t j_Ca40_K40_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ca40_K40_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ca40_K40_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ca40_K40_temp;

    AMREX_GPU_MANAGED table_t j_K40_Ar40_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_K40_Ar40_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_K40_Ar40_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_K40_Ar40_temp;

    AMREX_GPU_MANAGED table_t j_K40_Ca40_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_K40_Ca40_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_K40_Ca40_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_K40_Ca40_temp;

    AMREX_GPU_MANAGED table_t j_Ca41_K41_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ca41_K41_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ca41_K41_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ca41_K41_temp;

    AMREX_GPU_MANAGED table_t j_K41_Ca41_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_K41_Ca41_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_K41_Ca41_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_K41_Ca41_temp;

    AMREX_GPU_MANAGED table_t j_Ca43_Sc43_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ca43_Sc43_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ca43_Sc43_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ca43_Sc43_temp;

    AMREX_GPU_MANAGED table_t j_Sc43_Ca43_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Sc43_Ca43_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Sc43_Ca43_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Sc43_Ca43_temp;

    AMREX_GPU_MANAGED table_t j_Ca44_Sc44_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ca44_Sc44_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ca44_Sc44_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ca44_Sc44_temp;

    AMREX_GPU_MANAGED table_t j_Sc44_Ca44_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Sc44_Ca44_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Sc44_Ca44_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Sc44_Ca44_temp;

    AMREX_GPU_MANAGED table_t j_Sc44_Ti44_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Sc44_Ti44_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Sc44_Ti44_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Sc44_Ti44_temp;

    AMREX_GPU_MANAGED table_t j_Ti44_Sc44_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ti44_Sc44_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ti44_Sc44_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ti44_Sc44_temp;

    AMREX_GPU_MANAGED table_t j_Co53_Fe53_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Co53_Fe53_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Co53_Fe53_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Co53_Fe53_temp;

    AMREX_GPU_MANAGED table_t j_Fe53_Co53_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Fe53_Co53_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Fe53_Co53_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Fe53_Co53_temp;

    AMREX_GPU_MANAGED table_t j_Cu57_Ni57_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Cu57_Ni57_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Cu57_Ni57_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Cu57_Ni57_temp;

    AMREX_GPU_MANAGED table_t j_Ni57_Cu57_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ni57_Cu57_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ni57_Cu57_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ni57_Cu57_temp;

    AMREX_GPU_MANAGED table_t j_Ca45_Sc45_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ca45_Sc45_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ca45_Sc45_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ca45_Sc45_temp;

    AMREX_GPU_MANAGED table_t j_Sc45_Ca45_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Sc45_Ca45_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Sc45_Ca45_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Sc45_Ca45_temp;

    AMREX_GPU_MANAGED table_t j_Sc45_Ti45_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Sc45_Ti45_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Sc45_Ti45_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Sc45_Ti45_temp;

    AMREX_GPU_MANAGED table_t j_Ti45_Sc45_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ti45_Sc45_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ti45_Sc45_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ti45_Sc45_temp;

    AMREX_GPU_MANAGED table_t j_Ca46_Sc46_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ca46_Sc46_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ca46_Sc46_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ca46_Sc46_temp;

    AMREX_GPU_MANAGED table_t j_Sc46_Ca46_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Sc46_Ca46_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Sc46_Ca46_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Sc46_Ca46_temp;

    AMREX_GPU_MANAGED table_t j_Sc46_Ti46_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Sc46_Ti46_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Sc46_Ti46_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Sc46_Ti46_temp;

    AMREX_GPU_MANAGED table_t j_Ti46_Sc46_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ti46_Sc46_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ti46_Sc46_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ti46_Sc46_temp;

    AMREX_GPU_MANAGED table_t j_Ti46_V46_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ti46_V46_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ti46_V46_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ti46_V46_temp;

    AMREX_GPU_MANAGED table_t j_V46_Ti46_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_V46_Ti46_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_V46_Ti46_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_V46_Ti46_temp;

    AMREX_GPU_MANAGED table_t j_Ca47_Sc47_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ca47_Sc47_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ca47_Sc47_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ca47_Sc47_temp;

    AMREX_GPU_MANAGED table_t j_Sc47_Ca47_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Sc47_Ca47_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Sc47_Ca47_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Sc47_Ca47_temp;

    AMREX_GPU_MANAGED table_t j_Sc47_Ti47_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Sc47_Ti47_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Sc47_Ti47_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Sc47_Ti47_temp;

    AMREX_GPU_MANAGED table_t j_Ti47_Sc47_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ti47_Sc47_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ti47_Sc47_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ti47_Sc47_temp;

    AMREX_GPU_MANAGED table_t j_Ti47_V47_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ti47_V47_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ti47_V47_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ti47_V47_temp;

    AMREX_GPU_MANAGED table_t j_V47_Ti47_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_V47_Ti47_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_V47_Ti47_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_V47_Ti47_temp;

    AMREX_GPU_MANAGED table_t j_Ca48_Sc48_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ca48_Sc48_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ca48_Sc48_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ca48_Sc48_temp;

    AMREX_GPU_MANAGED table_t j_Cr48_V48_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Cr48_V48_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Cr48_V48_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Cr48_V48_temp;

    AMREX_GPU_MANAGED table_t j_Sc48_Ca48_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Sc48_Ca48_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Sc48_Ca48_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Sc48_Ca48_temp;

    AMREX_GPU_MANAGED table_t j_Sc48_Ti48_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Sc48_Ti48_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Sc48_Ti48_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Sc48_Ti48_temp;

    AMREX_GPU_MANAGED table_t j_Ti48_Sc48_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ti48_Sc48_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ti48_Sc48_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ti48_Sc48_temp;

    AMREX_GPU_MANAGED table_t j_Ti48_V48_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ti48_V48_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ti48_V48_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ti48_V48_temp;

    AMREX_GPU_MANAGED table_t j_V48_Cr48_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_V48_Cr48_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_V48_Cr48_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_V48_Cr48_temp;

    AMREX_GPU_MANAGED table_t j_V48_Ti48_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_V48_Ti48_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_V48_Ti48_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_V48_Ti48_temp;

    AMREX_GPU_MANAGED table_t j_Cr49_V49_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Cr49_V49_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Cr49_V49_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Cr49_V49_temp;

    AMREX_GPU_MANAGED table_t j_Sc49_Ti49_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Sc49_Ti49_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Sc49_Ti49_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Sc49_Ti49_temp;

    AMREX_GPU_MANAGED table_t j_Ti49_Sc49_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ti49_Sc49_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ti49_Sc49_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ti49_Sc49_temp;

    AMREX_GPU_MANAGED table_t j_Ti49_V49_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ti49_V49_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ti49_V49_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ti49_V49_temp;

    AMREX_GPU_MANAGED table_t j_V49_Cr49_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_V49_Cr49_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_V49_Cr49_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_V49_Cr49_temp;

    AMREX_GPU_MANAGED table_t j_V49_Ti49_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_V49_Ti49_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_V49_Ti49_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_V49_Ti49_temp;

    AMREX_GPU_MANAGED table_t j_Cr50_Mn50_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Cr50_Mn50_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Cr50_Mn50_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Cr50_Mn50_temp;

    AMREX_GPU_MANAGED table_t j_Cr50_V50_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Cr50_V50_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Cr50_V50_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Cr50_V50_temp;

    AMREX_GPU_MANAGED table_t j_Mn50_Cr50_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Mn50_Cr50_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Mn50_Cr50_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Mn50_Cr50_temp;

    AMREX_GPU_MANAGED table_t j_Ti50_V50_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ti50_V50_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ti50_V50_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ti50_V50_temp;

    AMREX_GPU_MANAGED table_t j_V50_Cr50_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_V50_Cr50_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_V50_Cr50_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_V50_Cr50_temp;

    AMREX_GPU_MANAGED table_t j_V50_Ti50_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_V50_Ti50_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_V50_Ti50_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_V50_Ti50_temp;

    AMREX_GPU_MANAGED table_t j_Cr51_Mn51_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Cr51_Mn51_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Cr51_Mn51_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Cr51_Mn51_temp;

    AMREX_GPU_MANAGED table_t j_Cr51_V51_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Cr51_V51_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Cr51_V51_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Cr51_V51_temp;

    AMREX_GPU_MANAGED table_t j_Mn51_Cr51_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Mn51_Cr51_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Mn51_Cr51_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Mn51_Cr51_temp;

    AMREX_GPU_MANAGED table_t j_Ti51_V51_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ti51_V51_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ti51_V51_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ti51_V51_temp;

    AMREX_GPU_MANAGED table_t j_V51_Cr51_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_V51_Cr51_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_V51_Cr51_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_V51_Cr51_temp;

    AMREX_GPU_MANAGED table_t j_V51_Ti51_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_V51_Ti51_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_V51_Ti51_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_V51_Ti51_temp;

    AMREX_GPU_MANAGED table_t j_Cr52_Mn52_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Cr52_Mn52_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Cr52_Mn52_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Cr52_Mn52_temp;

    AMREX_GPU_MANAGED table_t j_Cr52_V52_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Cr52_V52_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Cr52_V52_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Cr52_V52_temp;

    AMREX_GPU_MANAGED table_t j_Fe52_Mn52_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Fe52_Mn52_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Fe52_Mn52_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Fe52_Mn52_temp;

    AMREX_GPU_MANAGED table_t j_Mn52_Cr52_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Mn52_Cr52_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Mn52_Cr52_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Mn52_Cr52_temp;

    AMREX_GPU_MANAGED table_t j_Mn52_Fe52_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Mn52_Fe52_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Mn52_Fe52_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Mn52_Fe52_temp;

    AMREX_GPU_MANAGED table_t j_V52_Cr52_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_V52_Cr52_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_V52_Cr52_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_V52_Cr52_temp;

    AMREX_GPU_MANAGED table_t j_Cr53_Mn53_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Cr53_Mn53_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Cr53_Mn53_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Cr53_Mn53_temp;

    AMREX_GPU_MANAGED table_t j_Fe53_Mn53_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Fe53_Mn53_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Fe53_Mn53_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Fe53_Mn53_temp;

    AMREX_GPU_MANAGED table_t j_Mn53_Cr53_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Mn53_Cr53_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Mn53_Cr53_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Mn53_Cr53_temp;

    AMREX_GPU_MANAGED table_t j_Mn53_Fe53_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Mn53_Fe53_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Mn53_Fe53_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Mn53_Fe53_temp;

    AMREX_GPU_MANAGED table_t j_Co54_Fe54_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Co54_Fe54_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Co54_Fe54_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Co54_Fe54_temp;

    AMREX_GPU_MANAGED table_t j_Cr54_Mn54_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Cr54_Mn54_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Cr54_Mn54_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Cr54_Mn54_temp;

    AMREX_GPU_MANAGED table_t j_Fe54_Co54_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Fe54_Co54_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Fe54_Co54_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Fe54_Co54_temp;

    AMREX_GPU_MANAGED table_t j_Fe54_Mn54_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Fe54_Mn54_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Fe54_Mn54_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Fe54_Mn54_temp;

    AMREX_GPU_MANAGED table_t j_Mn54_Cr54_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Mn54_Cr54_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Mn54_Cr54_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Mn54_Cr54_temp;

    AMREX_GPU_MANAGED table_t j_Mn54_Fe54_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Mn54_Fe54_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Mn54_Fe54_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Mn54_Fe54_temp;

    AMREX_GPU_MANAGED table_t j_Co55_Fe55_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Co55_Fe55_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Co55_Fe55_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Co55_Fe55_temp;

    AMREX_GPU_MANAGED table_t j_Fe55_Co55_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Fe55_Co55_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Fe55_Co55_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Fe55_Co55_temp;

    AMREX_GPU_MANAGED table_t j_Fe55_Mn55_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Fe55_Mn55_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Fe55_Mn55_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Fe55_Mn55_temp;

    AMREX_GPU_MANAGED table_t j_Mn55_Fe55_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Mn55_Fe55_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Mn55_Fe55_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Mn55_Fe55_temp;

    AMREX_GPU_MANAGED table_t j_Co56_Fe56_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Co56_Fe56_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Co56_Fe56_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Co56_Fe56_temp;

    AMREX_GPU_MANAGED table_t j_Co56_Ni56_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Co56_Ni56_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Co56_Ni56_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Co56_Ni56_temp;

    AMREX_GPU_MANAGED table_t j_Fe56_Co56_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Fe56_Co56_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Fe56_Co56_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Fe56_Co56_temp;

    AMREX_GPU_MANAGED table_t j_Ni56_Co56_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ni56_Co56_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ni56_Co56_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ni56_Co56_temp;

    AMREX_GPU_MANAGED table_t j_Co57_Fe57_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Co57_Fe57_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Co57_Fe57_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Co57_Fe57_temp;

    AMREX_GPU_MANAGED table_t j_Co57_Ni57_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Co57_Ni57_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Co57_Ni57_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Co57_Ni57_temp;

    AMREX_GPU_MANAGED table_t j_Fe57_Co57_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Fe57_Co57_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Fe57_Co57_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Fe57_Co57_temp;

    AMREX_GPU_MANAGED table_t j_Ni57_Co57_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ni57_Co57_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ni57_Co57_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ni57_Co57_temp;

    AMREX_GPU_MANAGED table_t j_Co58_Fe58_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Co58_Fe58_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Co58_Fe58_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Co58_Fe58_temp;

    AMREX_GPU_MANAGED table_t j_Co58_Ni58_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Co58_Ni58_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Co58_Ni58_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Co58_Ni58_temp;

    AMREX_GPU_MANAGED table_t j_Cu58_Ni58_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Cu58_Ni58_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Cu58_Ni58_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Cu58_Ni58_temp;

    AMREX_GPU_MANAGED table_t j_Fe58_Co58_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Fe58_Co58_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Fe58_Co58_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Fe58_Co58_temp;

    AMREX_GPU_MANAGED table_t j_Ni58_Co58_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ni58_Co58_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ni58_Co58_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ni58_Co58_temp;

    AMREX_GPU_MANAGED table_t j_Ni58_Cu58_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ni58_Cu58_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ni58_Cu58_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ni58_Cu58_temp;

    AMREX_GPU_MANAGED table_t j_Co59_Ni59_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Co59_Ni59_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Co59_Ni59_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Co59_Ni59_temp;

    AMREX_GPU_MANAGED table_t j_Cu59_Ni59_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Cu59_Ni59_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Cu59_Ni59_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Cu59_Ni59_temp;

    AMREX_GPU_MANAGED table_t j_Ni59_Co59_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ni59_Co59_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ni59_Co59_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ni59_Co59_temp;

    AMREX_GPU_MANAGED table_t j_Ni59_Cu59_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ni59_Cu59_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ni59_Cu59_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ni59_Cu59_temp;

    AMREX_GPU_MANAGED table_t j_Cu60_Ni60_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Cu60_Ni60_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Cu60_Ni60_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Cu60_Ni60_temp;

    AMREX_GPU_MANAGED table_t j_Cu60_Zn60_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Cu60_Zn60_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Cu60_Zn60_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Cu60_Zn60_temp;

    AMREX_GPU_MANAGED table_t j_Ni60_Cu60_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ni60_Cu60_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ni60_Cu60_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ni60_Cu60_temp;

    AMREX_GPU_MANAGED table_t j_Zn60_Cu60_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Zn60_Cu60_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Zn60_Cu60_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Zn60_Cu60_temp;

    AMREX_GPU_MANAGED table_t j_Cu61_Ni61_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Cu61_Ni61_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Cu61_Ni61_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Cu61_Ni61_temp;

    AMREX_GPU_MANAGED table_t j_Cu61_Zn61_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Cu61_Zn61_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Cu61_Zn61_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Cu61_Zn61_temp;

    AMREX_GPU_MANAGED table_t j_Ni61_Cu61_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ni61_Cu61_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ni61_Cu61_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ni61_Cu61_temp;

    AMREX_GPU_MANAGED table_t j_Zn61_Cu61_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Zn61_Cu61_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Zn61_Cu61_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Zn61_Cu61_temp;

    AMREX_GPU_MANAGED table_t j_Cu62_Ni62_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Cu62_Ni62_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Cu62_Ni62_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Cu62_Ni62_temp;

    AMREX_GPU_MANAGED table_t j_Cu62_Zn62_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Cu62_Zn62_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Cu62_Zn62_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Cu62_Zn62_temp;

    AMREX_GPU_MANAGED table_t j_Ga62_Zn62_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ga62_Zn62_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ga62_Zn62_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ga62_Zn62_temp;

    AMREX_GPU_MANAGED table_t j_Ni62_Cu62_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ni62_Cu62_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ni62_Cu62_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ni62_Cu62_temp;

    AMREX_GPU_MANAGED table_t j_Zn62_Cu62_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Zn62_Cu62_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Zn62_Cu62_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Zn62_Cu62_temp;

    AMREX_GPU_MANAGED table_t j_Zn62_Ga62_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Zn62_Ga62_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Zn62_Ga62_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Zn62_Ga62_temp;

    AMREX_GPU_MANAGED table_t j_Cu63_Ni63_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Cu63_Ni63_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Cu63_Ni63_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Cu63_Ni63_temp;

    AMREX_GPU_MANAGED table_t j_Cu63_Zn63_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Cu63_Zn63_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Cu63_Zn63_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Cu63_Zn63_temp;

    AMREX_GPU_MANAGED table_t j_Ga63_Zn63_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ga63_Zn63_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ga63_Zn63_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ga63_Zn63_temp;

    AMREX_GPU_MANAGED table_t j_Ni63_Cu63_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ni63_Cu63_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ni63_Cu63_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ni63_Cu63_temp;

    AMREX_GPU_MANAGED table_t j_Zn63_Cu63_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Zn63_Cu63_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Zn63_Cu63_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Zn63_Cu63_temp;

    AMREX_GPU_MANAGED table_t j_Zn63_Ga63_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Zn63_Ga63_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Zn63_Ga63_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Zn63_Ga63_temp;

    AMREX_GPU_MANAGED table_t j_Cu64_Ni64_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Cu64_Ni64_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Cu64_Ni64_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Cu64_Ni64_temp;

    AMREX_GPU_MANAGED table_t j_Cu64_Zn64_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Cu64_Zn64_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Cu64_Zn64_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Cu64_Zn64_temp;

    AMREX_GPU_MANAGED table_t j_Ga64_Ge64_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ga64_Ge64_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ga64_Ge64_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ga64_Ge64_temp;

    AMREX_GPU_MANAGED table_t j_Ga64_Zn64_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ga64_Zn64_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ga64_Zn64_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ga64_Zn64_temp;

    AMREX_GPU_MANAGED table_t j_Ge64_Ga64_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ge64_Ga64_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ge64_Ga64_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ge64_Ga64_temp;

    AMREX_GPU_MANAGED table_t j_Ni64_Cu64_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ni64_Cu64_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ni64_Cu64_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ni64_Cu64_temp;

    AMREX_GPU_MANAGED table_t j_Zn64_Cu64_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Zn64_Cu64_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Zn64_Cu64_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Zn64_Cu64_temp;

    AMREX_GPU_MANAGED table_t j_Zn64_Ga64_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Zn64_Ga64_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Zn64_Ga64_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Zn64_Ga64_temp;

    AMREX_GPU_MANAGED table_t j_Cu65_Zn65_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Cu65_Zn65_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Cu65_Zn65_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Cu65_Zn65_temp;

    AMREX_GPU_MANAGED table_t j_Zn65_Cu65_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Zn65_Cu65_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Zn65_Cu65_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Zn65_Cu65_temp;

    AMREX_GPU_MANAGED table_t j_n_p_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_n_p_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_n_p_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_n_p_temp;

    AMREX_GPU_MANAGED table_t j_p_n_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_p_n_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_p_n_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_p_n_temp;


}


void init_tabular()
{

    amrex::Print() << "reading in network electron-capture / beta-decay tables..." << std::endl;

    using namespace rate_tables;

    j_Na21_Ne21_meta.ntemp = 13;
    j_Na21_Ne21_meta.nrhoy = 11;
    j_Na21_Ne21_meta.nvars = 6;
    j_Na21_Ne21_meta.nheader = 5;

    init_tab_info(j_Na21_Ne21_meta, "ffn-21na-21ne_electroncapture.dat", j_Na21_Ne21_rhoy, j_Na21_Ne21_temp, j_Na21_Ne21_data);


    j_Ne21_Na21_meta.ntemp = 13;
    j_Ne21_Na21_meta.nrhoy = 11;
    j_Ne21_Na21_meta.nvars = 6;
    j_Ne21_Na21_meta.nheader = 5;

    init_tab_info(j_Ne21_Na21_meta, "ffn-21ne-21na_betadecay.dat", j_Ne21_Na21_rhoy, j_Ne21_Na21_temp, j_Ne21_Na21_data);


    j_Na22_Ne22_meta.ntemp = 13;
    j_Na22_Ne22_meta.nrhoy = 11;
    j_Na22_Ne22_meta.nvars = 6;
    j_Na22_Ne22_meta.nheader = 5;

    init_tab_info(j_Na22_Ne22_meta, "ffn-22na-22ne_electroncapture.dat", j_Na22_Ne22_rhoy, j_Na22_Ne22_temp, j_Na22_Ne22_data);


    j_Ne22_Na22_meta.ntemp = 13;
    j_Ne22_Na22_meta.nrhoy = 11;
    j_Ne22_Na22_meta.nvars = 6;
    j_Ne22_Na22_meta.nheader = 5;

    init_tab_info(j_Ne22_Na22_meta, "ffn-22ne-22na_betadecay.dat", j_Ne22_Na22_rhoy, j_Ne22_Na22_temp, j_Ne22_Na22_data);


    j_Mg23_Na23_meta.ntemp = 13;
    j_Mg23_Na23_meta.nrhoy = 11;
    j_Mg23_Na23_meta.nvars = 6;
    j_Mg23_Na23_meta.nheader = 5;

    init_tab_info(j_Mg23_Na23_meta, "ffn-23mg-23na_electroncapture.dat", j_Mg23_Na23_rhoy, j_Mg23_Na23_temp, j_Mg23_Na23_data);


    j_Na23_Mg23_meta.ntemp = 13;
    j_Na23_Mg23_meta.nrhoy = 11;
    j_Na23_Mg23_meta.nvars = 6;
    j_Na23_Mg23_meta.nheader = 5;

    init_tab_info(j_Na23_Mg23_meta, "ffn-23na-23mg_betadecay.dat", j_Na23_Mg23_rhoy, j_Na23_Mg23_temp, j_Na23_Mg23_data);


    j_Al25_Mg25_meta.ntemp = 13;
    j_Al25_Mg25_meta.nrhoy = 11;
    j_Al25_Mg25_meta.nvars = 6;
    j_Al25_Mg25_meta.nheader = 5;

    init_tab_info(j_Al25_Mg25_meta, "ffn-25al-25mg_electroncapture.dat", j_Al25_Mg25_rhoy, j_Al25_Mg25_temp, j_Al25_Mg25_data);


    j_Mg25_Al25_meta.ntemp = 13;
    j_Mg25_Al25_meta.nrhoy = 11;
    j_Mg25_Al25_meta.nvars = 6;
    j_Mg25_Al25_meta.nheader = 5;

    init_tab_info(j_Mg25_Al25_meta, "ffn-25mg-25al_betadecay.dat", j_Mg25_Al25_rhoy, j_Mg25_Al25_temp, j_Mg25_Al25_data);


    j_Al26_Mg26_meta.ntemp = 13;
    j_Al26_Mg26_meta.nrhoy = 11;
    j_Al26_Mg26_meta.nvars = 6;
    j_Al26_Mg26_meta.nheader = 5;

    init_tab_info(j_Al26_Mg26_meta, "ffn-26al-26mg_electroncapture.dat", j_Al26_Mg26_rhoy, j_Al26_Mg26_temp, j_Al26_Mg26_data);


    j_Mg26_Al26_meta.ntemp = 13;
    j_Mg26_Al26_meta.nrhoy = 11;
    j_Mg26_Al26_meta.nvars = 6;
    j_Mg26_Al26_meta.nheader = 5;

    init_tab_info(j_Mg26_Al26_meta, "ffn-26mg-26al_betadecay.dat", j_Mg26_Al26_rhoy, j_Mg26_Al26_temp, j_Mg26_Al26_data);


    j_P29_Si29_meta.ntemp = 13;
    j_P29_Si29_meta.nrhoy = 11;
    j_P29_Si29_meta.nvars = 6;
    j_P29_Si29_meta.nheader = 5;

    init_tab_info(j_P29_Si29_meta, "ffn-29p-29si_electroncapture.dat", j_P29_Si29_rhoy, j_P29_Si29_temp, j_P29_Si29_data);


    j_Si29_P29_meta.ntemp = 13;
    j_Si29_P29_meta.nrhoy = 11;
    j_Si29_P29_meta.nvars = 6;
    j_Si29_P29_meta.nheader = 5;

    init_tab_info(j_Si29_P29_meta, "ffn-29si-29p_betadecay.dat", j_Si29_P29_rhoy, j_Si29_P29_temp, j_Si29_P29_data);


    j_P30_Si30_meta.ntemp = 13;
    j_P30_Si30_meta.nrhoy = 11;
    j_P30_Si30_meta.nvars = 6;
    j_P30_Si30_meta.nheader = 5;

    init_tab_info(j_P30_Si30_meta, "ffn-30p-30si_electroncapture.dat", j_P30_Si30_rhoy, j_P30_Si30_temp, j_P30_Si30_data);


    j_Si30_P30_meta.ntemp = 13;
    j_Si30_P30_meta.nrhoy = 11;
    j_Si30_P30_meta.nvars = 6;
    j_Si30_P30_meta.nheader = 5;

    init_tab_info(j_Si30_P30_meta, "ffn-30si-30p_betadecay.dat", j_Si30_P30_rhoy, j_Si30_P30_temp, j_Si30_P30_data);


    j_P31_Si31_meta.ntemp = 13;
    j_P31_Si31_meta.nrhoy = 11;
    j_P31_Si31_meta.nvars = 6;
    j_P31_Si31_meta.nheader = 5;

    init_tab_info(j_P31_Si31_meta, "ffn-31p-31si_electroncapture.dat", j_P31_Si31_rhoy, j_P31_Si31_temp, j_P31_Si31_data);


    j_Si31_P31_meta.ntemp = 13;
    j_Si31_P31_meta.nrhoy = 11;
    j_Si31_P31_meta.nvars = 6;
    j_Si31_P31_meta.nheader = 5;

    init_tab_info(j_Si31_P31_meta, "ffn-31si-31p_betadecay.dat", j_Si31_P31_rhoy, j_Si31_P31_temp, j_Si31_P31_data);


    j_P32_S32_meta.ntemp = 13;
    j_P32_S32_meta.nrhoy = 11;
    j_P32_S32_meta.nvars = 6;
    j_P32_S32_meta.nheader = 5;

    init_tab_info(j_P32_S32_meta, "ffn-32p-32s_betadecay.dat", j_P32_S32_rhoy, j_P32_S32_temp, j_P32_S32_data);


    j_P32_Si32_meta.ntemp = 13;
    j_P32_Si32_meta.nrhoy = 11;
    j_P32_Si32_meta.nvars = 6;
    j_P32_Si32_meta.nheader = 5;

    init_tab_info(j_P32_Si32_meta, "ffn-32p-32si_electroncapture.dat", j_P32_Si32_rhoy, j_P32_Si32_temp, j_P32_Si32_data);


    j_S32_P32_meta.ntemp = 13;
    j_S32_P32_meta.nrhoy = 11;
    j_S32_P32_meta.nvars = 6;
    j_S32_P32_meta.nheader = 5;

    init_tab_info(j_S32_P32_meta, "ffn-32s-32p_electroncapture.dat", j_S32_P32_rhoy, j_S32_P32_temp, j_S32_P32_data);


    j_Si32_P32_meta.ntemp = 13;
    j_Si32_P32_meta.nrhoy = 11;
    j_Si32_P32_meta.nvars = 6;
    j_Si32_P32_meta.nheader = 5;

    init_tab_info(j_Si32_P32_meta, "ffn-32si-32p_betadecay.dat", j_Si32_P32_rhoy, j_Si32_P32_temp, j_Si32_P32_data);


    j_Cl33_S33_meta.ntemp = 13;
    j_Cl33_S33_meta.nrhoy = 11;
    j_Cl33_S33_meta.nvars = 6;
    j_Cl33_S33_meta.nheader = 5;

    init_tab_info(j_Cl33_S33_meta, "ffn-33cl-33s_electroncapture.dat", j_Cl33_S33_rhoy, j_Cl33_S33_temp, j_Cl33_S33_data);


    j_P33_S33_meta.ntemp = 13;
    j_P33_S33_meta.nrhoy = 11;
    j_P33_S33_meta.nvars = 6;
    j_P33_S33_meta.nheader = 5;

    init_tab_info(j_P33_S33_meta, "ffn-33p-33s_betadecay.dat", j_P33_S33_rhoy, j_P33_S33_temp, j_P33_S33_data);


    j_S33_Cl33_meta.ntemp = 13;
    j_S33_Cl33_meta.nrhoy = 11;
    j_S33_Cl33_meta.nvars = 6;
    j_S33_Cl33_meta.nheader = 5;

    init_tab_info(j_S33_Cl33_meta, "ffn-33s-33cl_betadecay.dat", j_S33_Cl33_rhoy, j_S33_Cl33_temp, j_S33_Cl33_data);


    j_S33_P33_meta.ntemp = 13;
    j_S33_P33_meta.nrhoy = 11;
    j_S33_P33_meta.nvars = 6;
    j_S33_P33_meta.nheader = 5;

    init_tab_info(j_S33_P33_meta, "ffn-33s-33p_electroncapture.dat", j_S33_P33_rhoy, j_S33_P33_temp, j_S33_P33_data);


    j_Cl34_S34_meta.ntemp = 13;
    j_Cl34_S34_meta.nrhoy = 11;
    j_Cl34_S34_meta.nvars = 6;
    j_Cl34_S34_meta.nheader = 5;

    init_tab_info(j_Cl34_S34_meta, "ffn-34cl-34s_electroncapture.dat", j_Cl34_S34_rhoy, j_Cl34_S34_temp, j_Cl34_S34_data);


    j_S34_Cl34_meta.ntemp = 13;
    j_S34_Cl34_meta.nrhoy = 11;
    j_S34_Cl34_meta.nvars = 6;
    j_S34_Cl34_meta.nheader = 5;

    init_tab_info(j_S34_Cl34_meta, "ffn-34s-34cl_betadecay.dat", j_S34_Cl34_rhoy, j_S34_Cl34_temp, j_S34_Cl34_data);


    j_Cl35_S35_meta.ntemp = 13;
    j_Cl35_S35_meta.nrhoy = 11;
    j_Cl35_S35_meta.nvars = 6;
    j_Cl35_S35_meta.nheader = 5;

    init_tab_info(j_Cl35_S35_meta, "ffn-35cl-35s_electroncapture.dat", j_Cl35_S35_rhoy, j_Cl35_S35_temp, j_Cl35_S35_data);


    j_S35_Cl35_meta.ntemp = 13;
    j_S35_Cl35_meta.nrhoy = 11;
    j_S35_Cl35_meta.nvars = 6;
    j_S35_Cl35_meta.nheader = 5;

    init_tab_info(j_S35_Cl35_meta, "ffn-35s-35cl_betadecay.dat", j_S35_Cl35_rhoy, j_S35_Cl35_temp, j_S35_Cl35_data);


    j_Ar36_Cl36_meta.ntemp = 13;
    j_Ar36_Cl36_meta.nrhoy = 11;
    j_Ar36_Cl36_meta.nvars = 6;
    j_Ar36_Cl36_meta.nheader = 5;

    init_tab_info(j_Ar36_Cl36_meta, "ffn-36ar-36cl_electroncapture.dat", j_Ar36_Cl36_rhoy, j_Ar36_Cl36_temp, j_Ar36_Cl36_data);


    j_Cl36_Ar36_meta.ntemp = 13;
    j_Cl36_Ar36_meta.nrhoy = 11;
    j_Cl36_Ar36_meta.nvars = 6;
    j_Cl36_Ar36_meta.nheader = 5;

    init_tab_info(j_Cl36_Ar36_meta, "ffn-36cl-36ar_betadecay.dat", j_Cl36_Ar36_rhoy, j_Cl36_Ar36_temp, j_Cl36_Ar36_data);


    j_Cl36_S36_meta.ntemp = 13;
    j_Cl36_S36_meta.nrhoy = 11;
    j_Cl36_S36_meta.nvars = 6;
    j_Cl36_S36_meta.nheader = 5;

    init_tab_info(j_Cl36_S36_meta, "ffn-36cl-36s_electroncapture.dat", j_Cl36_S36_rhoy, j_Cl36_S36_temp, j_Cl36_S36_data);


    j_S36_Cl36_meta.ntemp = 13;
    j_S36_Cl36_meta.nrhoy = 11;
    j_S36_Cl36_meta.nvars = 6;
    j_S36_Cl36_meta.nheader = 5;

    init_tab_info(j_S36_Cl36_meta, "ffn-36s-36cl_betadecay.dat", j_S36_Cl36_rhoy, j_S36_Cl36_temp, j_S36_Cl36_data);


    j_Ar37_Cl37_meta.ntemp = 13;
    j_Ar37_Cl37_meta.nrhoy = 11;
    j_Ar37_Cl37_meta.nvars = 6;
    j_Ar37_Cl37_meta.nheader = 5;

    init_tab_info(j_Ar37_Cl37_meta, "ffn-37ar-37cl_electroncapture.dat", j_Ar37_Cl37_rhoy, j_Ar37_Cl37_temp, j_Ar37_Cl37_data);


    j_Ar37_K37_meta.ntemp = 13;
    j_Ar37_K37_meta.nrhoy = 11;
    j_Ar37_K37_meta.nvars = 6;
    j_Ar37_K37_meta.nheader = 5;

    init_tab_info(j_Ar37_K37_meta, "ffn-37ar-37k_betadecay.dat", j_Ar37_K37_rhoy, j_Ar37_K37_temp, j_Ar37_K37_data);


    j_Cl37_Ar37_meta.ntemp = 13;
    j_Cl37_Ar37_meta.nrhoy = 11;
    j_Cl37_Ar37_meta.nvars = 6;
    j_Cl37_Ar37_meta.nheader = 5;

    init_tab_info(j_Cl37_Ar37_meta, "ffn-37cl-37ar_betadecay.dat", j_Cl37_Ar37_rhoy, j_Cl37_Ar37_temp, j_Cl37_Ar37_data);


    j_K37_Ar37_meta.ntemp = 13;
    j_K37_Ar37_meta.nrhoy = 11;
    j_K37_Ar37_meta.nvars = 6;
    j_K37_Ar37_meta.nheader = 5;

    init_tab_info(j_K37_Ar37_meta, "ffn-37k-37ar_electroncapture.dat", j_K37_Ar37_rhoy, j_K37_Ar37_temp, j_K37_Ar37_data);


    j_Ar38_K38_meta.ntemp = 13;
    j_Ar38_K38_meta.nrhoy = 11;
    j_Ar38_K38_meta.nvars = 6;
    j_Ar38_K38_meta.nheader = 5;

    init_tab_info(j_Ar38_K38_meta, "ffn-38ar-38k_betadecay.dat", j_Ar38_K38_rhoy, j_Ar38_K38_temp, j_Ar38_K38_data);


    j_K38_Ar38_meta.ntemp = 13;
    j_K38_Ar38_meta.nrhoy = 11;
    j_K38_Ar38_meta.nvars = 6;
    j_K38_Ar38_meta.nheader = 5;

    init_tab_info(j_K38_Ar38_meta, "ffn-38k-38ar_electroncapture.dat", j_K38_Ar38_rhoy, j_K38_Ar38_temp, j_K38_Ar38_data);


    j_Ar39_K39_meta.ntemp = 13;
    j_Ar39_K39_meta.nrhoy = 11;
    j_Ar39_K39_meta.nvars = 6;
    j_Ar39_K39_meta.nheader = 5;

    init_tab_info(j_Ar39_K39_meta, "ffn-39ar-39k_betadecay.dat", j_Ar39_K39_rhoy, j_Ar39_K39_temp, j_Ar39_K39_data);


    j_K39_Ar39_meta.ntemp = 13;
    j_K39_Ar39_meta.nrhoy = 11;
    j_K39_Ar39_meta.nvars = 6;
    j_K39_Ar39_meta.nheader = 5;

    init_tab_info(j_K39_Ar39_meta, "ffn-39k-39ar_electroncapture.dat", j_K39_Ar39_rhoy, j_K39_Ar39_temp, j_K39_Ar39_data);


    j_Ar40_K40_meta.ntemp = 13;
    j_Ar40_K40_meta.nrhoy = 11;
    j_Ar40_K40_meta.nvars = 6;
    j_Ar40_K40_meta.nheader = 5;

    init_tab_info(j_Ar40_K40_meta, "ffn-40ar-40k_betadecay.dat", j_Ar40_K40_rhoy, j_Ar40_K40_temp, j_Ar40_K40_data);


    j_Ca40_K40_meta.ntemp = 13;
    j_Ca40_K40_meta.nrhoy = 11;
    j_Ca40_K40_meta.nvars = 6;
    j_Ca40_K40_meta.nheader = 5;

    init_tab_info(j_Ca40_K40_meta, "ffn-40ca-40k_electroncapture.dat", j_Ca40_K40_rhoy, j_Ca40_K40_temp, j_Ca40_K40_data);


    j_K40_Ar40_meta.ntemp = 13;
    j_K40_Ar40_meta.nrhoy = 11;
    j_K40_Ar40_meta.nvars = 6;
    j_K40_Ar40_meta.nheader = 5;

    init_tab_info(j_K40_Ar40_meta, "ffn-40k-40ar_electroncapture.dat", j_K40_Ar40_rhoy, j_K40_Ar40_temp, j_K40_Ar40_data);


    j_K40_Ca40_meta.ntemp = 13;
    j_K40_Ca40_meta.nrhoy = 11;
    j_K40_Ca40_meta.nvars = 6;
    j_K40_Ca40_meta.nheader = 5;

    init_tab_info(j_K40_Ca40_meta, "ffn-40k-40ca_betadecay.dat", j_K40_Ca40_rhoy, j_K40_Ca40_temp, j_K40_Ca40_data);


    j_Ca41_K41_meta.ntemp = 13;
    j_Ca41_K41_meta.nrhoy = 11;
    j_Ca41_K41_meta.nvars = 6;
    j_Ca41_K41_meta.nheader = 5;

    init_tab_info(j_Ca41_K41_meta, "ffn-41ca-41k_electroncapture.dat", j_Ca41_K41_rhoy, j_Ca41_K41_temp, j_Ca41_K41_data);


    j_K41_Ca41_meta.ntemp = 13;
    j_K41_Ca41_meta.nrhoy = 11;
    j_K41_Ca41_meta.nvars = 6;
    j_K41_Ca41_meta.nheader = 5;

    init_tab_info(j_K41_Ca41_meta, "ffn-41k-41ca_betadecay.dat", j_K41_Ca41_rhoy, j_K41_Ca41_temp, j_K41_Ca41_data);


    j_Ca43_Sc43_meta.ntemp = 13;
    j_Ca43_Sc43_meta.nrhoy = 11;
    j_Ca43_Sc43_meta.nvars = 6;
    j_Ca43_Sc43_meta.nheader = 5;

    init_tab_info(j_Ca43_Sc43_meta, "ffn-43ca-43sc_betadecay.dat", j_Ca43_Sc43_rhoy, j_Ca43_Sc43_temp, j_Ca43_Sc43_data);


    j_Sc43_Ca43_meta.ntemp = 13;
    j_Sc43_Ca43_meta.nrhoy = 11;
    j_Sc43_Ca43_meta.nvars = 6;
    j_Sc43_Ca43_meta.nheader = 5;

    init_tab_info(j_Sc43_Ca43_meta, "ffn-43sc-43ca_electroncapture.dat", j_Sc43_Ca43_rhoy, j_Sc43_Ca43_temp, j_Sc43_Ca43_data);


    j_Ca44_Sc44_meta.ntemp = 13;
    j_Ca44_Sc44_meta.nrhoy = 11;
    j_Ca44_Sc44_meta.nvars = 6;
    j_Ca44_Sc44_meta.nheader = 5;

    init_tab_info(j_Ca44_Sc44_meta, "ffn-44ca-44sc_betadecay.dat", j_Ca44_Sc44_rhoy, j_Ca44_Sc44_temp, j_Ca44_Sc44_data);


    j_Sc44_Ca44_meta.ntemp = 13;
    j_Sc44_Ca44_meta.nrhoy = 11;
    j_Sc44_Ca44_meta.nvars = 6;
    j_Sc44_Ca44_meta.nheader = 5;

    init_tab_info(j_Sc44_Ca44_meta, "ffn-44sc-44ca_electroncapture.dat", j_Sc44_Ca44_rhoy, j_Sc44_Ca44_temp, j_Sc44_Ca44_data);


    j_Sc44_Ti44_meta.ntemp = 13;
    j_Sc44_Ti44_meta.nrhoy = 11;
    j_Sc44_Ti44_meta.nvars = 6;
    j_Sc44_Ti44_meta.nheader = 5;

    init_tab_info(j_Sc44_Ti44_meta, "ffn-44sc-44ti_betadecay.dat", j_Sc44_Ti44_rhoy, j_Sc44_Ti44_temp, j_Sc44_Ti44_data);


    j_Ti44_Sc44_meta.ntemp = 13;
    j_Ti44_Sc44_meta.nrhoy = 11;
    j_Ti44_Sc44_meta.nvars = 6;
    j_Ti44_Sc44_meta.nheader = 5;

    init_tab_info(j_Ti44_Sc44_meta, "ffn-44ti-44sc_electroncapture.dat", j_Ti44_Sc44_rhoy, j_Ti44_Sc44_temp, j_Ti44_Sc44_data);


    j_Co53_Fe53_meta.ntemp = 13;
    j_Co53_Fe53_meta.nrhoy = 11;
    j_Co53_Fe53_meta.nvars = 6;
    j_Co53_Fe53_meta.nheader = 5;

    init_tab_info(j_Co53_Fe53_meta, "ffn-53co-53fe_electroncapture.dat", j_Co53_Fe53_rhoy, j_Co53_Fe53_temp, j_Co53_Fe53_data);


    j_Fe53_Co53_meta.ntemp = 13;
    j_Fe53_Co53_meta.nrhoy = 11;
    j_Fe53_Co53_meta.nvars = 6;
    j_Fe53_Co53_meta.nheader = 5;

    init_tab_info(j_Fe53_Co53_meta, "ffn-53fe-53co_betadecay.dat", j_Fe53_Co53_rhoy, j_Fe53_Co53_temp, j_Fe53_Co53_data);


    j_Cu57_Ni57_meta.ntemp = 13;
    j_Cu57_Ni57_meta.nrhoy = 11;
    j_Cu57_Ni57_meta.nvars = 6;
    j_Cu57_Ni57_meta.nheader = 5;

    init_tab_info(j_Cu57_Ni57_meta, "ffn-57cu-57ni_electroncapture.dat", j_Cu57_Ni57_rhoy, j_Cu57_Ni57_temp, j_Cu57_Ni57_data);


    j_Ni57_Cu57_meta.ntemp = 13;
    j_Ni57_Cu57_meta.nrhoy = 11;
    j_Ni57_Cu57_meta.nvars = 6;
    j_Ni57_Cu57_meta.nheader = 5;

    init_tab_info(j_Ni57_Cu57_meta, "ffn-57ni-57cu_betadecay.dat", j_Ni57_Cu57_rhoy, j_Ni57_Cu57_temp, j_Ni57_Cu57_data);


    j_Ca45_Sc45_meta.ntemp = 13;
    j_Ca45_Sc45_meta.nrhoy = 11;
    j_Ca45_Sc45_meta.nvars = 6;
    j_Ca45_Sc45_meta.nheader = 5;

    init_tab_info(j_Ca45_Sc45_meta, "langanke-45ca-45sc_betadecay.dat", j_Ca45_Sc45_rhoy, j_Ca45_Sc45_temp, j_Ca45_Sc45_data);


    j_Sc45_Ca45_meta.ntemp = 13;
    j_Sc45_Ca45_meta.nrhoy = 11;
    j_Sc45_Ca45_meta.nvars = 6;
    j_Sc45_Ca45_meta.nheader = 5;

    init_tab_info(j_Sc45_Ca45_meta, "langanke-45sc-45ca_electroncapture.dat", j_Sc45_Ca45_rhoy, j_Sc45_Ca45_temp, j_Sc45_Ca45_data);


    j_Sc45_Ti45_meta.ntemp = 13;
    j_Sc45_Ti45_meta.nrhoy = 11;
    j_Sc45_Ti45_meta.nvars = 6;
    j_Sc45_Ti45_meta.nheader = 5;

    init_tab_info(j_Sc45_Ti45_meta, "langanke-45sc-45ti_betadecay.dat", j_Sc45_Ti45_rhoy, j_Sc45_Ti45_temp, j_Sc45_Ti45_data);


    j_Ti45_Sc45_meta.ntemp = 13;
    j_Ti45_Sc45_meta.nrhoy = 11;
    j_Ti45_Sc45_meta.nvars = 6;
    j_Ti45_Sc45_meta.nheader = 5;

    init_tab_info(j_Ti45_Sc45_meta, "langanke-45ti-45sc_electroncapture.dat", j_Ti45_Sc45_rhoy, j_Ti45_Sc45_temp, j_Ti45_Sc45_data);


    j_Ca46_Sc46_meta.ntemp = 13;
    j_Ca46_Sc46_meta.nrhoy = 11;
    j_Ca46_Sc46_meta.nvars = 6;
    j_Ca46_Sc46_meta.nheader = 5;

    init_tab_info(j_Ca46_Sc46_meta, "langanke-46ca-46sc_betadecay.dat", j_Ca46_Sc46_rhoy, j_Ca46_Sc46_temp, j_Ca46_Sc46_data);


    j_Sc46_Ca46_meta.ntemp = 13;
    j_Sc46_Ca46_meta.nrhoy = 11;
    j_Sc46_Ca46_meta.nvars = 6;
    j_Sc46_Ca46_meta.nheader = 5;

    init_tab_info(j_Sc46_Ca46_meta, "langanke-46sc-46ca_electroncapture.dat", j_Sc46_Ca46_rhoy, j_Sc46_Ca46_temp, j_Sc46_Ca46_data);


    j_Sc46_Ti46_meta.ntemp = 13;
    j_Sc46_Ti46_meta.nrhoy = 11;
    j_Sc46_Ti46_meta.nvars = 6;
    j_Sc46_Ti46_meta.nheader = 5;

    init_tab_info(j_Sc46_Ti46_meta, "langanke-46sc-46ti_betadecay.dat", j_Sc46_Ti46_rhoy, j_Sc46_Ti46_temp, j_Sc46_Ti46_data);


    j_Ti46_Sc46_meta.ntemp = 13;
    j_Ti46_Sc46_meta.nrhoy = 11;
    j_Ti46_Sc46_meta.nvars = 6;
    j_Ti46_Sc46_meta.nheader = 5;

    init_tab_info(j_Ti46_Sc46_meta, "langanke-46ti-46sc_electroncapture.dat", j_Ti46_Sc46_rhoy, j_Ti46_Sc46_temp, j_Ti46_Sc46_data);


    j_Ti46_V46_meta.ntemp = 13;
    j_Ti46_V46_meta.nrhoy = 11;
    j_Ti46_V46_meta.nvars = 6;
    j_Ti46_V46_meta.nheader = 5;

    init_tab_info(j_Ti46_V46_meta, "langanke-46ti-46v_betadecay.dat", j_Ti46_V46_rhoy, j_Ti46_V46_temp, j_Ti46_V46_data);


    j_V46_Ti46_meta.ntemp = 13;
    j_V46_Ti46_meta.nrhoy = 11;
    j_V46_Ti46_meta.nvars = 6;
    j_V46_Ti46_meta.nheader = 5;

    init_tab_info(j_V46_Ti46_meta, "langanke-46v-46ti_electroncapture.dat", j_V46_Ti46_rhoy, j_V46_Ti46_temp, j_V46_Ti46_data);


    j_Ca47_Sc47_meta.ntemp = 13;
    j_Ca47_Sc47_meta.nrhoy = 11;
    j_Ca47_Sc47_meta.nvars = 6;
    j_Ca47_Sc47_meta.nheader = 5;

    init_tab_info(j_Ca47_Sc47_meta, "langanke-47ca-47sc_betadecay.dat", j_Ca47_Sc47_rhoy, j_Ca47_Sc47_temp, j_Ca47_Sc47_data);


    j_Sc47_Ca47_meta.ntemp = 13;
    j_Sc47_Ca47_meta.nrhoy = 11;
    j_Sc47_Ca47_meta.nvars = 6;
    j_Sc47_Ca47_meta.nheader = 5;

    init_tab_info(j_Sc47_Ca47_meta, "langanke-47sc-47ca_electroncapture.dat", j_Sc47_Ca47_rhoy, j_Sc47_Ca47_temp, j_Sc47_Ca47_data);


    j_Sc47_Ti47_meta.ntemp = 13;
    j_Sc47_Ti47_meta.nrhoy = 11;
    j_Sc47_Ti47_meta.nvars = 6;
    j_Sc47_Ti47_meta.nheader = 5;

    init_tab_info(j_Sc47_Ti47_meta, "langanke-47sc-47ti_betadecay.dat", j_Sc47_Ti47_rhoy, j_Sc47_Ti47_temp, j_Sc47_Ti47_data);


    j_Ti47_Sc47_meta.ntemp = 13;
    j_Ti47_Sc47_meta.nrhoy = 11;
    j_Ti47_Sc47_meta.nvars = 6;
    j_Ti47_Sc47_meta.nheader = 5;

    init_tab_info(j_Ti47_Sc47_meta, "langanke-47ti-47sc_electroncapture.dat", j_Ti47_Sc47_rhoy, j_Ti47_Sc47_temp, j_Ti47_Sc47_data);


    j_Ti47_V47_meta.ntemp = 13;
    j_Ti47_V47_meta.nrhoy = 11;
    j_Ti47_V47_meta.nvars = 6;
    j_Ti47_V47_meta.nheader = 5;

    init_tab_info(j_Ti47_V47_meta, "langanke-47ti-47v_betadecay.dat", j_Ti47_V47_rhoy, j_Ti47_V47_temp, j_Ti47_V47_data);


    j_V47_Ti47_meta.ntemp = 13;
    j_V47_Ti47_meta.nrhoy = 11;
    j_V47_Ti47_meta.nvars = 6;
    j_V47_Ti47_meta.nheader = 5;

    init_tab_info(j_V47_Ti47_meta, "langanke-47v-47ti_electroncapture.dat", j_V47_Ti47_rhoy, j_V47_Ti47_temp, j_V47_Ti47_data);


    j_Ca48_Sc48_meta.ntemp = 13;
    j_Ca48_Sc48_meta.nrhoy = 11;
    j_Ca48_Sc48_meta.nvars = 6;
    j_Ca48_Sc48_meta.nheader = 5;

    init_tab_info(j_Ca48_Sc48_meta, "langanke-48ca-48sc_betadecay.dat", j_Ca48_Sc48_rhoy, j_Ca48_Sc48_temp, j_Ca48_Sc48_data);


    j_Cr48_V48_meta.ntemp = 13;
    j_Cr48_V48_meta.nrhoy = 11;
    j_Cr48_V48_meta.nvars = 6;
    j_Cr48_V48_meta.nheader = 5;

    init_tab_info(j_Cr48_V48_meta, "langanke-48cr-48v_electroncapture.dat", j_Cr48_V48_rhoy, j_Cr48_V48_temp, j_Cr48_V48_data);


    j_Sc48_Ca48_meta.ntemp = 13;
    j_Sc48_Ca48_meta.nrhoy = 11;
    j_Sc48_Ca48_meta.nvars = 6;
    j_Sc48_Ca48_meta.nheader = 5;

    init_tab_info(j_Sc48_Ca48_meta, "langanke-48sc-48ca_electroncapture.dat", j_Sc48_Ca48_rhoy, j_Sc48_Ca48_temp, j_Sc48_Ca48_data);


    j_Sc48_Ti48_meta.ntemp = 13;
    j_Sc48_Ti48_meta.nrhoy = 11;
    j_Sc48_Ti48_meta.nvars = 6;
    j_Sc48_Ti48_meta.nheader = 5;

    init_tab_info(j_Sc48_Ti48_meta, "langanke-48sc-48ti_betadecay.dat", j_Sc48_Ti48_rhoy, j_Sc48_Ti48_temp, j_Sc48_Ti48_data);


    j_Ti48_Sc48_meta.ntemp = 13;
    j_Ti48_Sc48_meta.nrhoy = 11;
    j_Ti48_Sc48_meta.nvars = 6;
    j_Ti48_Sc48_meta.nheader = 5;

    init_tab_info(j_Ti48_Sc48_meta, "langanke-48ti-48sc_electroncapture.dat", j_Ti48_Sc48_rhoy, j_Ti48_Sc48_temp, j_Ti48_Sc48_data);


    j_Ti48_V48_meta.ntemp = 13;
    j_Ti48_V48_meta.nrhoy = 11;
    j_Ti48_V48_meta.nvars = 6;
    j_Ti48_V48_meta.nheader = 5;

    init_tab_info(j_Ti48_V48_meta, "langanke-48ti-48v_betadecay.dat", j_Ti48_V48_rhoy, j_Ti48_V48_temp, j_Ti48_V48_data);


    j_V48_Cr48_meta.ntemp = 13;
    j_V48_Cr48_meta.nrhoy = 11;
    j_V48_Cr48_meta.nvars = 6;
    j_V48_Cr48_meta.nheader = 5;

    init_tab_info(j_V48_Cr48_meta, "langanke-48v-48cr_betadecay.dat", j_V48_Cr48_rhoy, j_V48_Cr48_temp, j_V48_Cr48_data);


    j_V48_Ti48_meta.ntemp = 13;
    j_V48_Ti48_meta.nrhoy = 11;
    j_V48_Ti48_meta.nvars = 6;
    j_V48_Ti48_meta.nheader = 5;

    init_tab_info(j_V48_Ti48_meta, "langanke-48v-48ti_electroncapture.dat", j_V48_Ti48_rhoy, j_V48_Ti48_temp, j_V48_Ti48_data);


    j_Cr49_V49_meta.ntemp = 13;
    j_Cr49_V49_meta.nrhoy = 11;
    j_Cr49_V49_meta.nvars = 6;
    j_Cr49_V49_meta.nheader = 5;

    init_tab_info(j_Cr49_V49_meta, "langanke-49cr-49v_electroncapture.dat", j_Cr49_V49_rhoy, j_Cr49_V49_temp, j_Cr49_V49_data);


    j_Sc49_Ti49_meta.ntemp = 13;
    j_Sc49_Ti49_meta.nrhoy = 11;
    j_Sc49_Ti49_meta.nvars = 6;
    j_Sc49_Ti49_meta.nheader = 5;

    init_tab_info(j_Sc49_Ti49_meta, "langanke-49sc-49ti_betadecay.dat", j_Sc49_Ti49_rhoy, j_Sc49_Ti49_temp, j_Sc49_Ti49_data);


    j_Ti49_Sc49_meta.ntemp = 13;
    j_Ti49_Sc49_meta.nrhoy = 11;
    j_Ti49_Sc49_meta.nvars = 6;
    j_Ti49_Sc49_meta.nheader = 5;

    init_tab_info(j_Ti49_Sc49_meta, "langanke-49ti-49sc_electroncapture.dat", j_Ti49_Sc49_rhoy, j_Ti49_Sc49_temp, j_Ti49_Sc49_data);


    j_Ti49_V49_meta.ntemp = 13;
    j_Ti49_V49_meta.nrhoy = 11;
    j_Ti49_V49_meta.nvars = 6;
    j_Ti49_V49_meta.nheader = 5;

    init_tab_info(j_Ti49_V49_meta, "langanke-49ti-49v_betadecay.dat", j_Ti49_V49_rhoy, j_Ti49_V49_temp, j_Ti49_V49_data);


    j_V49_Cr49_meta.ntemp = 13;
    j_V49_Cr49_meta.nrhoy = 11;
    j_V49_Cr49_meta.nvars = 6;
    j_V49_Cr49_meta.nheader = 5;

    init_tab_info(j_V49_Cr49_meta, "langanke-49v-49cr_betadecay.dat", j_V49_Cr49_rhoy, j_V49_Cr49_temp, j_V49_Cr49_data);


    j_V49_Ti49_meta.ntemp = 13;
    j_V49_Ti49_meta.nrhoy = 11;
    j_V49_Ti49_meta.nvars = 6;
    j_V49_Ti49_meta.nheader = 5;

    init_tab_info(j_V49_Ti49_meta, "langanke-49v-49ti_electroncapture.dat", j_V49_Ti49_rhoy, j_V49_Ti49_temp, j_V49_Ti49_data);


    j_Cr50_Mn50_meta.ntemp = 13;
    j_Cr50_Mn50_meta.nrhoy = 11;
    j_Cr50_Mn50_meta.nvars = 6;
    j_Cr50_Mn50_meta.nheader = 5;

    init_tab_info(j_Cr50_Mn50_meta, "langanke-50cr-50mn_betadecay.dat", j_Cr50_Mn50_rhoy, j_Cr50_Mn50_temp, j_Cr50_Mn50_data);


    j_Cr50_V50_meta.ntemp = 13;
    j_Cr50_V50_meta.nrhoy = 11;
    j_Cr50_V50_meta.nvars = 6;
    j_Cr50_V50_meta.nheader = 5;

    init_tab_info(j_Cr50_V50_meta, "langanke-50cr-50v_electroncapture.dat", j_Cr50_V50_rhoy, j_Cr50_V50_temp, j_Cr50_V50_data);


    j_Mn50_Cr50_meta.ntemp = 13;
    j_Mn50_Cr50_meta.nrhoy = 11;
    j_Mn50_Cr50_meta.nvars = 6;
    j_Mn50_Cr50_meta.nheader = 5;

    init_tab_info(j_Mn50_Cr50_meta, "langanke-50mn-50cr_electroncapture.dat", j_Mn50_Cr50_rhoy, j_Mn50_Cr50_temp, j_Mn50_Cr50_data);


    j_Ti50_V50_meta.ntemp = 13;
    j_Ti50_V50_meta.nrhoy = 11;
    j_Ti50_V50_meta.nvars = 6;
    j_Ti50_V50_meta.nheader = 5;

    init_tab_info(j_Ti50_V50_meta, "langanke-50ti-50v_betadecay.dat", j_Ti50_V50_rhoy, j_Ti50_V50_temp, j_Ti50_V50_data);


    j_V50_Cr50_meta.ntemp = 13;
    j_V50_Cr50_meta.nrhoy = 11;
    j_V50_Cr50_meta.nvars = 6;
    j_V50_Cr50_meta.nheader = 5;

    init_tab_info(j_V50_Cr50_meta, "langanke-50v-50cr_betadecay.dat", j_V50_Cr50_rhoy, j_V50_Cr50_temp, j_V50_Cr50_data);


    j_V50_Ti50_meta.ntemp = 13;
    j_V50_Ti50_meta.nrhoy = 11;
    j_V50_Ti50_meta.nvars = 6;
    j_V50_Ti50_meta.nheader = 5;

    init_tab_info(j_V50_Ti50_meta, "langanke-50v-50ti_electroncapture.dat", j_V50_Ti50_rhoy, j_V50_Ti50_temp, j_V50_Ti50_data);


    j_Cr51_Mn51_meta.ntemp = 13;
    j_Cr51_Mn51_meta.nrhoy = 11;
    j_Cr51_Mn51_meta.nvars = 6;
    j_Cr51_Mn51_meta.nheader = 5;

    init_tab_info(j_Cr51_Mn51_meta, "langanke-51cr-51mn_betadecay.dat", j_Cr51_Mn51_rhoy, j_Cr51_Mn51_temp, j_Cr51_Mn51_data);


    j_Cr51_V51_meta.ntemp = 13;
    j_Cr51_V51_meta.nrhoy = 11;
    j_Cr51_V51_meta.nvars = 6;
    j_Cr51_V51_meta.nheader = 5;

    init_tab_info(j_Cr51_V51_meta, "langanke-51cr-51v_electroncapture.dat", j_Cr51_V51_rhoy, j_Cr51_V51_temp, j_Cr51_V51_data);


    j_Mn51_Cr51_meta.ntemp = 13;
    j_Mn51_Cr51_meta.nrhoy = 11;
    j_Mn51_Cr51_meta.nvars = 6;
    j_Mn51_Cr51_meta.nheader = 5;

    init_tab_info(j_Mn51_Cr51_meta, "langanke-51mn-51cr_electroncapture.dat", j_Mn51_Cr51_rhoy, j_Mn51_Cr51_temp, j_Mn51_Cr51_data);


    j_Ti51_V51_meta.ntemp = 13;
    j_Ti51_V51_meta.nrhoy = 11;
    j_Ti51_V51_meta.nvars = 6;
    j_Ti51_V51_meta.nheader = 5;

    init_tab_info(j_Ti51_V51_meta, "langanke-51ti-51v_betadecay.dat", j_Ti51_V51_rhoy, j_Ti51_V51_temp, j_Ti51_V51_data);


    j_V51_Cr51_meta.ntemp = 13;
    j_V51_Cr51_meta.nrhoy = 11;
    j_V51_Cr51_meta.nvars = 6;
    j_V51_Cr51_meta.nheader = 5;

    init_tab_info(j_V51_Cr51_meta, "langanke-51v-51cr_betadecay.dat", j_V51_Cr51_rhoy, j_V51_Cr51_temp, j_V51_Cr51_data);


    j_V51_Ti51_meta.ntemp = 13;
    j_V51_Ti51_meta.nrhoy = 11;
    j_V51_Ti51_meta.nvars = 6;
    j_V51_Ti51_meta.nheader = 5;

    init_tab_info(j_V51_Ti51_meta, "langanke-51v-51ti_electroncapture.dat", j_V51_Ti51_rhoy, j_V51_Ti51_temp, j_V51_Ti51_data);


    j_Cr52_Mn52_meta.ntemp = 13;
    j_Cr52_Mn52_meta.nrhoy = 11;
    j_Cr52_Mn52_meta.nvars = 6;
    j_Cr52_Mn52_meta.nheader = 5;

    init_tab_info(j_Cr52_Mn52_meta, "langanke-52cr-52mn_betadecay.dat", j_Cr52_Mn52_rhoy, j_Cr52_Mn52_temp, j_Cr52_Mn52_data);


    j_Cr52_V52_meta.ntemp = 13;
    j_Cr52_V52_meta.nrhoy = 11;
    j_Cr52_V52_meta.nvars = 6;
    j_Cr52_V52_meta.nheader = 5;

    init_tab_info(j_Cr52_V52_meta, "langanke-52cr-52v_electroncapture.dat", j_Cr52_V52_rhoy, j_Cr52_V52_temp, j_Cr52_V52_data);


    j_Fe52_Mn52_meta.ntemp = 13;
    j_Fe52_Mn52_meta.nrhoy = 11;
    j_Fe52_Mn52_meta.nvars = 6;
    j_Fe52_Mn52_meta.nheader = 5;

    init_tab_info(j_Fe52_Mn52_meta, "langanke-52fe-52mn_electroncapture.dat", j_Fe52_Mn52_rhoy, j_Fe52_Mn52_temp, j_Fe52_Mn52_data);


    j_Mn52_Cr52_meta.ntemp = 13;
    j_Mn52_Cr52_meta.nrhoy = 11;
    j_Mn52_Cr52_meta.nvars = 6;
    j_Mn52_Cr52_meta.nheader = 5;

    init_tab_info(j_Mn52_Cr52_meta, "langanke-52mn-52cr_electroncapture.dat", j_Mn52_Cr52_rhoy, j_Mn52_Cr52_temp, j_Mn52_Cr52_data);


    j_Mn52_Fe52_meta.ntemp = 13;
    j_Mn52_Fe52_meta.nrhoy = 11;
    j_Mn52_Fe52_meta.nvars = 6;
    j_Mn52_Fe52_meta.nheader = 5;

    init_tab_info(j_Mn52_Fe52_meta, "langanke-52mn-52fe_betadecay.dat", j_Mn52_Fe52_rhoy, j_Mn52_Fe52_temp, j_Mn52_Fe52_data);


    j_V52_Cr52_meta.ntemp = 13;
    j_V52_Cr52_meta.nrhoy = 11;
    j_V52_Cr52_meta.nvars = 6;
    j_V52_Cr52_meta.nheader = 5;

    init_tab_info(j_V52_Cr52_meta, "langanke-52v-52cr_betadecay.dat", j_V52_Cr52_rhoy, j_V52_Cr52_temp, j_V52_Cr52_data);


    j_Cr53_Mn53_meta.ntemp = 13;
    j_Cr53_Mn53_meta.nrhoy = 11;
    j_Cr53_Mn53_meta.nvars = 6;
    j_Cr53_Mn53_meta.nheader = 5;

    init_tab_info(j_Cr53_Mn53_meta, "langanke-53cr-53mn_betadecay.dat", j_Cr53_Mn53_rhoy, j_Cr53_Mn53_temp, j_Cr53_Mn53_data);


    j_Fe53_Mn53_meta.ntemp = 13;
    j_Fe53_Mn53_meta.nrhoy = 11;
    j_Fe53_Mn53_meta.nvars = 6;
    j_Fe53_Mn53_meta.nheader = 5;

    init_tab_info(j_Fe53_Mn53_meta, "langanke-53fe-53mn_electroncapture.dat", j_Fe53_Mn53_rhoy, j_Fe53_Mn53_temp, j_Fe53_Mn53_data);


    j_Mn53_Cr53_meta.ntemp = 13;
    j_Mn53_Cr53_meta.nrhoy = 11;
    j_Mn53_Cr53_meta.nvars = 6;
    j_Mn53_Cr53_meta.nheader = 5;

    init_tab_info(j_Mn53_Cr53_meta, "langanke-53mn-53cr_electroncapture.dat", j_Mn53_Cr53_rhoy, j_Mn53_Cr53_temp, j_Mn53_Cr53_data);


    j_Mn53_Fe53_meta.ntemp = 13;
    j_Mn53_Fe53_meta.nrhoy = 11;
    j_Mn53_Fe53_meta.nvars = 6;
    j_Mn53_Fe53_meta.nheader = 5;

    init_tab_info(j_Mn53_Fe53_meta, "langanke-53mn-53fe_betadecay.dat", j_Mn53_Fe53_rhoy, j_Mn53_Fe53_temp, j_Mn53_Fe53_data);


    j_Co54_Fe54_meta.ntemp = 13;
    j_Co54_Fe54_meta.nrhoy = 11;
    j_Co54_Fe54_meta.nvars = 6;
    j_Co54_Fe54_meta.nheader = 5;

    init_tab_info(j_Co54_Fe54_meta, "langanke-54co-54fe_electroncapture.dat", j_Co54_Fe54_rhoy, j_Co54_Fe54_temp, j_Co54_Fe54_data);


    j_Cr54_Mn54_meta.ntemp = 13;
    j_Cr54_Mn54_meta.nrhoy = 11;
    j_Cr54_Mn54_meta.nvars = 6;
    j_Cr54_Mn54_meta.nheader = 5;

    init_tab_info(j_Cr54_Mn54_meta, "langanke-54cr-54mn_betadecay.dat", j_Cr54_Mn54_rhoy, j_Cr54_Mn54_temp, j_Cr54_Mn54_data);


    j_Fe54_Co54_meta.ntemp = 13;
    j_Fe54_Co54_meta.nrhoy = 11;
    j_Fe54_Co54_meta.nvars = 6;
    j_Fe54_Co54_meta.nheader = 5;

    init_tab_info(j_Fe54_Co54_meta, "langanke-54fe-54co_betadecay.dat", j_Fe54_Co54_rhoy, j_Fe54_Co54_temp, j_Fe54_Co54_data);


    j_Fe54_Mn54_meta.ntemp = 13;
    j_Fe54_Mn54_meta.nrhoy = 11;
    j_Fe54_Mn54_meta.nvars = 6;
    j_Fe54_Mn54_meta.nheader = 5;

    init_tab_info(j_Fe54_Mn54_meta, "langanke-54fe-54mn_electroncapture.dat", j_Fe54_Mn54_rhoy, j_Fe54_Mn54_temp, j_Fe54_Mn54_data);


    j_Mn54_Cr54_meta.ntemp = 13;
    j_Mn54_Cr54_meta.nrhoy = 11;
    j_Mn54_Cr54_meta.nvars = 6;
    j_Mn54_Cr54_meta.nheader = 5;

    init_tab_info(j_Mn54_Cr54_meta, "langanke-54mn-54cr_electroncapture.dat", j_Mn54_Cr54_rhoy, j_Mn54_Cr54_temp, j_Mn54_Cr54_data);


    j_Mn54_Fe54_meta.ntemp = 13;
    j_Mn54_Fe54_meta.nrhoy = 11;
    j_Mn54_Fe54_meta.nvars = 6;
    j_Mn54_Fe54_meta.nheader = 5;

    init_tab_info(j_Mn54_Fe54_meta, "langanke-54mn-54fe_betadecay.dat", j_Mn54_Fe54_rhoy, j_Mn54_Fe54_temp, j_Mn54_Fe54_data);


    j_Co55_Fe55_meta.ntemp = 13;
    j_Co55_Fe55_meta.nrhoy = 11;
    j_Co55_Fe55_meta.nvars = 6;
    j_Co55_Fe55_meta.nheader = 5;

    init_tab_info(j_Co55_Fe55_meta, "langanke-55co-55fe_electroncapture.dat", j_Co55_Fe55_rhoy, j_Co55_Fe55_temp, j_Co55_Fe55_data);


    j_Fe55_Co55_meta.ntemp = 13;
    j_Fe55_Co55_meta.nrhoy = 11;
    j_Fe55_Co55_meta.nvars = 6;
    j_Fe55_Co55_meta.nheader = 5;

    init_tab_info(j_Fe55_Co55_meta, "langanke-55fe-55co_betadecay.dat", j_Fe55_Co55_rhoy, j_Fe55_Co55_temp, j_Fe55_Co55_data);


    j_Fe55_Mn55_meta.ntemp = 13;
    j_Fe55_Mn55_meta.nrhoy = 11;
    j_Fe55_Mn55_meta.nvars = 6;
    j_Fe55_Mn55_meta.nheader = 5;

    init_tab_info(j_Fe55_Mn55_meta, "langanke-55fe-55mn_electroncapture.dat", j_Fe55_Mn55_rhoy, j_Fe55_Mn55_temp, j_Fe55_Mn55_data);


    j_Mn55_Fe55_meta.ntemp = 13;
    j_Mn55_Fe55_meta.nrhoy = 11;
    j_Mn55_Fe55_meta.nvars = 6;
    j_Mn55_Fe55_meta.nheader = 5;

    init_tab_info(j_Mn55_Fe55_meta, "langanke-55mn-55fe_betadecay.dat", j_Mn55_Fe55_rhoy, j_Mn55_Fe55_temp, j_Mn55_Fe55_data);


    j_Co56_Fe56_meta.ntemp = 13;
    j_Co56_Fe56_meta.nrhoy = 11;
    j_Co56_Fe56_meta.nvars = 6;
    j_Co56_Fe56_meta.nheader = 5;

    init_tab_info(j_Co56_Fe56_meta, "langanke-56co-56fe_electroncapture.dat", j_Co56_Fe56_rhoy, j_Co56_Fe56_temp, j_Co56_Fe56_data);


    j_Co56_Ni56_meta.ntemp = 13;
    j_Co56_Ni56_meta.nrhoy = 11;
    j_Co56_Ni56_meta.nvars = 6;
    j_Co56_Ni56_meta.nheader = 5;

    init_tab_info(j_Co56_Ni56_meta, "langanke-56co-56ni_betadecay.dat", j_Co56_Ni56_rhoy, j_Co56_Ni56_temp, j_Co56_Ni56_data);


    j_Fe56_Co56_meta.ntemp = 13;
    j_Fe56_Co56_meta.nrhoy = 11;
    j_Fe56_Co56_meta.nvars = 6;
    j_Fe56_Co56_meta.nheader = 5;

    init_tab_info(j_Fe56_Co56_meta, "langanke-56fe-56co_betadecay.dat", j_Fe56_Co56_rhoy, j_Fe56_Co56_temp, j_Fe56_Co56_data);


    j_Ni56_Co56_meta.ntemp = 13;
    j_Ni56_Co56_meta.nrhoy = 11;
    j_Ni56_Co56_meta.nvars = 6;
    j_Ni56_Co56_meta.nheader = 5;

    init_tab_info(j_Ni56_Co56_meta, "langanke-56ni-56co_electroncapture.dat", j_Ni56_Co56_rhoy, j_Ni56_Co56_temp, j_Ni56_Co56_data);


    j_Co57_Fe57_meta.ntemp = 13;
    j_Co57_Fe57_meta.nrhoy = 11;
    j_Co57_Fe57_meta.nvars = 6;
    j_Co57_Fe57_meta.nheader = 5;

    init_tab_info(j_Co57_Fe57_meta, "langanke-57co-57fe_electroncapture.dat", j_Co57_Fe57_rhoy, j_Co57_Fe57_temp, j_Co57_Fe57_data);


    j_Co57_Ni57_meta.ntemp = 13;
    j_Co57_Ni57_meta.nrhoy = 11;
    j_Co57_Ni57_meta.nvars = 6;
    j_Co57_Ni57_meta.nheader = 5;

    init_tab_info(j_Co57_Ni57_meta, "langanke-57co-57ni_betadecay.dat", j_Co57_Ni57_rhoy, j_Co57_Ni57_temp, j_Co57_Ni57_data);


    j_Fe57_Co57_meta.ntemp = 13;
    j_Fe57_Co57_meta.nrhoy = 11;
    j_Fe57_Co57_meta.nvars = 6;
    j_Fe57_Co57_meta.nheader = 5;

    init_tab_info(j_Fe57_Co57_meta, "langanke-57fe-57co_betadecay.dat", j_Fe57_Co57_rhoy, j_Fe57_Co57_temp, j_Fe57_Co57_data);


    j_Ni57_Co57_meta.ntemp = 13;
    j_Ni57_Co57_meta.nrhoy = 11;
    j_Ni57_Co57_meta.nvars = 6;
    j_Ni57_Co57_meta.nheader = 5;

    init_tab_info(j_Ni57_Co57_meta, "langanke-57ni-57co_electroncapture.dat", j_Ni57_Co57_rhoy, j_Ni57_Co57_temp, j_Ni57_Co57_data);


    j_Co58_Fe58_meta.ntemp = 13;
    j_Co58_Fe58_meta.nrhoy = 11;
    j_Co58_Fe58_meta.nvars = 6;
    j_Co58_Fe58_meta.nheader = 5;

    init_tab_info(j_Co58_Fe58_meta, "langanke-58co-58fe_electroncapture.dat", j_Co58_Fe58_rhoy, j_Co58_Fe58_temp, j_Co58_Fe58_data);


    j_Co58_Ni58_meta.ntemp = 13;
    j_Co58_Ni58_meta.nrhoy = 11;
    j_Co58_Ni58_meta.nvars = 6;
    j_Co58_Ni58_meta.nheader = 5;

    init_tab_info(j_Co58_Ni58_meta, "langanke-58co-58ni_betadecay.dat", j_Co58_Ni58_rhoy, j_Co58_Ni58_temp, j_Co58_Ni58_data);


    j_Cu58_Ni58_meta.ntemp = 13;
    j_Cu58_Ni58_meta.nrhoy = 11;
    j_Cu58_Ni58_meta.nvars = 6;
    j_Cu58_Ni58_meta.nheader = 5;

    init_tab_info(j_Cu58_Ni58_meta, "langanke-58cu-58ni_electroncapture.dat", j_Cu58_Ni58_rhoy, j_Cu58_Ni58_temp, j_Cu58_Ni58_data);


    j_Fe58_Co58_meta.ntemp = 13;
    j_Fe58_Co58_meta.nrhoy = 11;
    j_Fe58_Co58_meta.nvars = 6;
    j_Fe58_Co58_meta.nheader = 5;

    init_tab_info(j_Fe58_Co58_meta, "langanke-58fe-58co_betadecay.dat", j_Fe58_Co58_rhoy, j_Fe58_Co58_temp, j_Fe58_Co58_data);


    j_Ni58_Co58_meta.ntemp = 13;
    j_Ni58_Co58_meta.nrhoy = 11;
    j_Ni58_Co58_meta.nvars = 6;
    j_Ni58_Co58_meta.nheader = 5;

    init_tab_info(j_Ni58_Co58_meta, "langanke-58ni-58co_electroncapture.dat", j_Ni58_Co58_rhoy, j_Ni58_Co58_temp, j_Ni58_Co58_data);


    j_Ni58_Cu58_meta.ntemp = 13;
    j_Ni58_Cu58_meta.nrhoy = 11;
    j_Ni58_Cu58_meta.nvars = 6;
    j_Ni58_Cu58_meta.nheader = 5;

    init_tab_info(j_Ni58_Cu58_meta, "langanke-58ni-58cu_betadecay.dat", j_Ni58_Cu58_rhoy, j_Ni58_Cu58_temp, j_Ni58_Cu58_data);


    j_Co59_Ni59_meta.ntemp = 13;
    j_Co59_Ni59_meta.nrhoy = 11;
    j_Co59_Ni59_meta.nvars = 6;
    j_Co59_Ni59_meta.nheader = 5;

    init_tab_info(j_Co59_Ni59_meta, "langanke-59co-59ni_betadecay.dat", j_Co59_Ni59_rhoy, j_Co59_Ni59_temp, j_Co59_Ni59_data);


    j_Cu59_Ni59_meta.ntemp = 13;
    j_Cu59_Ni59_meta.nrhoy = 11;
    j_Cu59_Ni59_meta.nvars = 6;
    j_Cu59_Ni59_meta.nheader = 5;

    init_tab_info(j_Cu59_Ni59_meta, "langanke-59cu-59ni_electroncapture.dat", j_Cu59_Ni59_rhoy, j_Cu59_Ni59_temp, j_Cu59_Ni59_data);


    j_Ni59_Co59_meta.ntemp = 13;
    j_Ni59_Co59_meta.nrhoy = 11;
    j_Ni59_Co59_meta.nvars = 6;
    j_Ni59_Co59_meta.nheader = 5;

    init_tab_info(j_Ni59_Co59_meta, "langanke-59ni-59co_electroncapture.dat", j_Ni59_Co59_rhoy, j_Ni59_Co59_temp, j_Ni59_Co59_data);


    j_Ni59_Cu59_meta.ntemp = 13;
    j_Ni59_Cu59_meta.nrhoy = 11;
    j_Ni59_Cu59_meta.nvars = 6;
    j_Ni59_Cu59_meta.nheader = 5;

    init_tab_info(j_Ni59_Cu59_meta, "langanke-59ni-59cu_betadecay.dat", j_Ni59_Cu59_rhoy, j_Ni59_Cu59_temp, j_Ni59_Cu59_data);


    j_Cu60_Ni60_meta.ntemp = 13;
    j_Cu60_Ni60_meta.nrhoy = 11;
    j_Cu60_Ni60_meta.nvars = 6;
    j_Cu60_Ni60_meta.nheader = 5;

    init_tab_info(j_Cu60_Ni60_meta, "langanke-60cu-60ni_electroncapture.dat", j_Cu60_Ni60_rhoy, j_Cu60_Ni60_temp, j_Cu60_Ni60_data);


    j_Cu60_Zn60_meta.ntemp = 13;
    j_Cu60_Zn60_meta.nrhoy = 11;
    j_Cu60_Zn60_meta.nvars = 6;
    j_Cu60_Zn60_meta.nheader = 5;

    init_tab_info(j_Cu60_Zn60_meta, "langanke-60cu-60zn_betadecay.dat", j_Cu60_Zn60_rhoy, j_Cu60_Zn60_temp, j_Cu60_Zn60_data);


    j_Ni60_Cu60_meta.ntemp = 13;
    j_Ni60_Cu60_meta.nrhoy = 11;
    j_Ni60_Cu60_meta.nvars = 6;
    j_Ni60_Cu60_meta.nheader = 5;

    init_tab_info(j_Ni60_Cu60_meta, "langanke-60ni-60cu_betadecay.dat", j_Ni60_Cu60_rhoy, j_Ni60_Cu60_temp, j_Ni60_Cu60_data);


    j_Zn60_Cu60_meta.ntemp = 13;
    j_Zn60_Cu60_meta.nrhoy = 11;
    j_Zn60_Cu60_meta.nvars = 6;
    j_Zn60_Cu60_meta.nheader = 5;

    init_tab_info(j_Zn60_Cu60_meta, "langanke-60zn-60cu_electroncapture.dat", j_Zn60_Cu60_rhoy, j_Zn60_Cu60_temp, j_Zn60_Cu60_data);


    j_Cu61_Ni61_meta.ntemp = 13;
    j_Cu61_Ni61_meta.nrhoy = 11;
    j_Cu61_Ni61_meta.nvars = 6;
    j_Cu61_Ni61_meta.nheader = 5;

    init_tab_info(j_Cu61_Ni61_meta, "langanke-61cu-61ni_electroncapture.dat", j_Cu61_Ni61_rhoy, j_Cu61_Ni61_temp, j_Cu61_Ni61_data);


    j_Cu61_Zn61_meta.ntemp = 13;
    j_Cu61_Zn61_meta.nrhoy = 11;
    j_Cu61_Zn61_meta.nvars = 6;
    j_Cu61_Zn61_meta.nheader = 5;

    init_tab_info(j_Cu61_Zn61_meta, "langanke-61cu-61zn_betadecay.dat", j_Cu61_Zn61_rhoy, j_Cu61_Zn61_temp, j_Cu61_Zn61_data);


    j_Ni61_Cu61_meta.ntemp = 13;
    j_Ni61_Cu61_meta.nrhoy = 11;
    j_Ni61_Cu61_meta.nvars = 6;
    j_Ni61_Cu61_meta.nheader = 5;

    init_tab_info(j_Ni61_Cu61_meta, "langanke-61ni-61cu_betadecay.dat", j_Ni61_Cu61_rhoy, j_Ni61_Cu61_temp, j_Ni61_Cu61_data);


    j_Zn61_Cu61_meta.ntemp = 13;
    j_Zn61_Cu61_meta.nrhoy = 11;
    j_Zn61_Cu61_meta.nvars = 6;
    j_Zn61_Cu61_meta.nheader = 5;

    init_tab_info(j_Zn61_Cu61_meta, "langanke-61zn-61cu_electroncapture.dat", j_Zn61_Cu61_rhoy, j_Zn61_Cu61_temp, j_Zn61_Cu61_data);


    j_Cu62_Ni62_meta.ntemp = 13;
    j_Cu62_Ni62_meta.nrhoy = 11;
    j_Cu62_Ni62_meta.nvars = 6;
    j_Cu62_Ni62_meta.nheader = 5;

    init_tab_info(j_Cu62_Ni62_meta, "langanke-62cu-62ni_electroncapture.dat", j_Cu62_Ni62_rhoy, j_Cu62_Ni62_temp, j_Cu62_Ni62_data);


    j_Cu62_Zn62_meta.ntemp = 13;
    j_Cu62_Zn62_meta.nrhoy = 11;
    j_Cu62_Zn62_meta.nvars = 6;
    j_Cu62_Zn62_meta.nheader = 5;

    init_tab_info(j_Cu62_Zn62_meta, "langanke-62cu-62zn_betadecay.dat", j_Cu62_Zn62_rhoy, j_Cu62_Zn62_temp, j_Cu62_Zn62_data);


    j_Ga62_Zn62_meta.ntemp = 13;
    j_Ga62_Zn62_meta.nrhoy = 11;
    j_Ga62_Zn62_meta.nvars = 6;
    j_Ga62_Zn62_meta.nheader = 5;

    init_tab_info(j_Ga62_Zn62_meta, "langanke-62ga-62zn_electroncapture.dat", j_Ga62_Zn62_rhoy, j_Ga62_Zn62_temp, j_Ga62_Zn62_data);


    j_Ni62_Cu62_meta.ntemp = 13;
    j_Ni62_Cu62_meta.nrhoy = 11;
    j_Ni62_Cu62_meta.nvars = 6;
    j_Ni62_Cu62_meta.nheader = 5;

    init_tab_info(j_Ni62_Cu62_meta, "langanke-62ni-62cu_betadecay.dat", j_Ni62_Cu62_rhoy, j_Ni62_Cu62_temp, j_Ni62_Cu62_data);


    j_Zn62_Cu62_meta.ntemp = 13;
    j_Zn62_Cu62_meta.nrhoy = 11;
    j_Zn62_Cu62_meta.nvars = 6;
    j_Zn62_Cu62_meta.nheader = 5;

    init_tab_info(j_Zn62_Cu62_meta, "langanke-62zn-62cu_electroncapture.dat", j_Zn62_Cu62_rhoy, j_Zn62_Cu62_temp, j_Zn62_Cu62_data);


    j_Zn62_Ga62_meta.ntemp = 13;
    j_Zn62_Ga62_meta.nrhoy = 11;
    j_Zn62_Ga62_meta.nvars = 6;
    j_Zn62_Ga62_meta.nheader = 5;

    init_tab_info(j_Zn62_Ga62_meta, "langanke-62zn-62ga_betadecay.dat", j_Zn62_Ga62_rhoy, j_Zn62_Ga62_temp, j_Zn62_Ga62_data);


    j_Cu63_Ni63_meta.ntemp = 13;
    j_Cu63_Ni63_meta.nrhoy = 11;
    j_Cu63_Ni63_meta.nvars = 6;
    j_Cu63_Ni63_meta.nheader = 5;

    init_tab_info(j_Cu63_Ni63_meta, "langanke-63cu-63ni_electroncapture.dat", j_Cu63_Ni63_rhoy, j_Cu63_Ni63_temp, j_Cu63_Ni63_data);


    j_Cu63_Zn63_meta.ntemp = 13;
    j_Cu63_Zn63_meta.nrhoy = 11;
    j_Cu63_Zn63_meta.nvars = 6;
    j_Cu63_Zn63_meta.nheader = 5;

    init_tab_info(j_Cu63_Zn63_meta, "langanke-63cu-63zn_betadecay.dat", j_Cu63_Zn63_rhoy, j_Cu63_Zn63_temp, j_Cu63_Zn63_data);


    j_Ga63_Zn63_meta.ntemp = 13;
    j_Ga63_Zn63_meta.nrhoy = 11;
    j_Ga63_Zn63_meta.nvars = 6;
    j_Ga63_Zn63_meta.nheader = 5;

    init_tab_info(j_Ga63_Zn63_meta, "langanke-63ga-63zn_electroncapture.dat", j_Ga63_Zn63_rhoy, j_Ga63_Zn63_temp, j_Ga63_Zn63_data);


    j_Ni63_Cu63_meta.ntemp = 13;
    j_Ni63_Cu63_meta.nrhoy = 11;
    j_Ni63_Cu63_meta.nvars = 6;
    j_Ni63_Cu63_meta.nheader = 5;

    init_tab_info(j_Ni63_Cu63_meta, "langanke-63ni-63cu_betadecay.dat", j_Ni63_Cu63_rhoy, j_Ni63_Cu63_temp, j_Ni63_Cu63_data);


    j_Zn63_Cu63_meta.ntemp = 13;
    j_Zn63_Cu63_meta.nrhoy = 11;
    j_Zn63_Cu63_meta.nvars = 6;
    j_Zn63_Cu63_meta.nheader = 5;

    init_tab_info(j_Zn63_Cu63_meta, "langanke-63zn-63cu_electroncapture.dat", j_Zn63_Cu63_rhoy, j_Zn63_Cu63_temp, j_Zn63_Cu63_data);


    j_Zn63_Ga63_meta.ntemp = 13;
    j_Zn63_Ga63_meta.nrhoy = 11;
    j_Zn63_Ga63_meta.nvars = 6;
    j_Zn63_Ga63_meta.nheader = 5;

    init_tab_info(j_Zn63_Ga63_meta, "langanke-63zn-63ga_betadecay.dat", j_Zn63_Ga63_rhoy, j_Zn63_Ga63_temp, j_Zn63_Ga63_data);


    j_Cu64_Ni64_meta.ntemp = 13;
    j_Cu64_Ni64_meta.nrhoy = 11;
    j_Cu64_Ni64_meta.nvars = 6;
    j_Cu64_Ni64_meta.nheader = 5;

    init_tab_info(j_Cu64_Ni64_meta, "langanke-64cu-64ni_electroncapture.dat", j_Cu64_Ni64_rhoy, j_Cu64_Ni64_temp, j_Cu64_Ni64_data);


    j_Cu64_Zn64_meta.ntemp = 13;
    j_Cu64_Zn64_meta.nrhoy = 11;
    j_Cu64_Zn64_meta.nvars = 6;
    j_Cu64_Zn64_meta.nheader = 5;

    init_tab_info(j_Cu64_Zn64_meta, "langanke-64cu-64zn_betadecay.dat", j_Cu64_Zn64_rhoy, j_Cu64_Zn64_temp, j_Cu64_Zn64_data);


    j_Ga64_Ge64_meta.ntemp = 13;
    j_Ga64_Ge64_meta.nrhoy = 11;
    j_Ga64_Ge64_meta.nvars = 6;
    j_Ga64_Ge64_meta.nheader = 5;

    init_tab_info(j_Ga64_Ge64_meta, "langanke-64ga-64ge_betadecay.dat", j_Ga64_Ge64_rhoy, j_Ga64_Ge64_temp, j_Ga64_Ge64_data);


    j_Ga64_Zn64_meta.ntemp = 13;
    j_Ga64_Zn64_meta.nrhoy = 11;
    j_Ga64_Zn64_meta.nvars = 6;
    j_Ga64_Zn64_meta.nheader = 5;

    init_tab_info(j_Ga64_Zn64_meta, "langanke-64ga-64zn_electroncapture.dat", j_Ga64_Zn64_rhoy, j_Ga64_Zn64_temp, j_Ga64_Zn64_data);


    j_Ge64_Ga64_meta.ntemp = 13;
    j_Ge64_Ga64_meta.nrhoy = 11;
    j_Ge64_Ga64_meta.nvars = 6;
    j_Ge64_Ga64_meta.nheader = 5;

    init_tab_info(j_Ge64_Ga64_meta, "langanke-64ge-64ga_electroncapture.dat", j_Ge64_Ga64_rhoy, j_Ge64_Ga64_temp, j_Ge64_Ga64_data);


    j_Ni64_Cu64_meta.ntemp = 13;
    j_Ni64_Cu64_meta.nrhoy = 11;
    j_Ni64_Cu64_meta.nvars = 6;
    j_Ni64_Cu64_meta.nheader = 5;

    init_tab_info(j_Ni64_Cu64_meta, "langanke-64ni-64cu_betadecay.dat", j_Ni64_Cu64_rhoy, j_Ni64_Cu64_temp, j_Ni64_Cu64_data);


    j_Zn64_Cu64_meta.ntemp = 13;
    j_Zn64_Cu64_meta.nrhoy = 11;
    j_Zn64_Cu64_meta.nvars = 6;
    j_Zn64_Cu64_meta.nheader = 5;

    init_tab_info(j_Zn64_Cu64_meta, "langanke-64zn-64cu_electroncapture.dat", j_Zn64_Cu64_rhoy, j_Zn64_Cu64_temp, j_Zn64_Cu64_data);


    j_Zn64_Ga64_meta.ntemp = 13;
    j_Zn64_Ga64_meta.nrhoy = 11;
    j_Zn64_Ga64_meta.nvars = 6;
    j_Zn64_Ga64_meta.nheader = 5;

    init_tab_info(j_Zn64_Ga64_meta, "langanke-64zn-64ga_betadecay.dat", j_Zn64_Ga64_rhoy, j_Zn64_Ga64_temp, j_Zn64_Ga64_data);


    j_Cu65_Zn65_meta.ntemp = 13;
    j_Cu65_Zn65_meta.nrhoy = 11;
    j_Cu65_Zn65_meta.nvars = 6;
    j_Cu65_Zn65_meta.nheader = 5;

    init_tab_info(j_Cu65_Zn65_meta, "langanke-65cu-65zn_betadecay.dat", j_Cu65_Zn65_rhoy, j_Cu65_Zn65_temp, j_Cu65_Zn65_data);


    j_Zn65_Cu65_meta.ntemp = 13;
    j_Zn65_Cu65_meta.nrhoy = 11;
    j_Zn65_Cu65_meta.nvars = 6;
    j_Zn65_Cu65_meta.nheader = 5;

    init_tab_info(j_Zn65_Cu65_meta, "langanke-65zn-65cu_electroncapture.dat", j_Zn65_Cu65_rhoy, j_Zn65_Cu65_temp, j_Zn65_Cu65_data);


    j_n_p_meta.ntemp = 13;
    j_n_p_meta.nrhoy = 11;
    j_n_p_meta.nvars = 6;
    j_n_p_meta.nheader = 5;

    init_tab_info(j_n_p_meta, "langanke-n-p_betadecay.dat", j_n_p_rhoy, j_n_p_temp, j_n_p_data);


    j_p_n_meta.ntemp = 13;
    j_p_n_meta.nrhoy = 11;
    j_p_n_meta.nvars = 6;
    j_p_n_meta.nheader = 5;

    init_tab_info(j_p_n_meta, "langanke-p-n_electroncapture.dat", j_p_n_rhoy, j_p_n_temp, j_p_n_data);



}
