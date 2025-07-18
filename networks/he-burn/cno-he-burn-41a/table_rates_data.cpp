#include <AMReX_Array.H>
#include <string>
#include <table_rates.H>
#include <AMReX_Print.H>

namespace rate_tables
{

    AMREX_GPU_MANAGED table_t j_F17_O17_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 12, 1, 11, 1, 6> j_F17_O17_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_F17_O17_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 12> j_F17_O17_temp;

    AMREX_GPU_MANAGED table_t j_O17_F17_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 12, 1, 11, 1, 6> j_O17_F17_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_O17_F17_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 12> j_O17_F17_temp;

    AMREX_GPU_MANAGED table_t j_F18_Ne18_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 12, 1, 11, 1, 6> j_F18_Ne18_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_F18_Ne18_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 12> j_F18_Ne18_temp;

    AMREX_GPU_MANAGED table_t j_F18_O18_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 12, 1, 11, 1, 6> j_F18_O18_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_F18_O18_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 12> j_F18_O18_temp;

    AMREX_GPU_MANAGED table_t j_Ne18_F18_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 12, 1, 11, 1, 6> j_Ne18_F18_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ne18_F18_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 12> j_Ne18_F18_temp;

    AMREX_GPU_MANAGED table_t j_O18_F18_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 12, 1, 11, 1, 6> j_O18_F18_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_O18_F18_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 12> j_O18_F18_temp;

    AMREX_GPU_MANAGED table_t j_F19_Ne19_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 12, 1, 11, 1, 6> j_F19_Ne19_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_F19_Ne19_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 12> j_F19_Ne19_temp;

    AMREX_GPU_MANAGED table_t j_Ne19_F19_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 12, 1, 11, 1, 6> j_Ne19_F19_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ne19_F19_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 12> j_Ne19_F19_temp;

    AMREX_GPU_MANAGED table_t j_Na21_Ne21_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 12, 1, 11, 1, 6> j_Na21_Ne21_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Na21_Ne21_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 12> j_Na21_Ne21_temp;

    AMREX_GPU_MANAGED table_t j_Ne21_Na21_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 12, 1, 11, 1, 6> j_Ne21_Na21_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ne21_Na21_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 12> j_Ne21_Na21_temp;

    AMREX_GPU_MANAGED table_t j_Mg22_Na22_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 12, 1, 11, 1, 6> j_Mg22_Na22_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Mg22_Na22_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 12> j_Mg22_Na22_temp;

    AMREX_GPU_MANAGED table_t j_Na22_Mg22_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 12, 1, 11, 1, 6> j_Na22_Mg22_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Na22_Mg22_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 12> j_Na22_Mg22_temp;

    AMREX_GPU_MANAGED table_t j_Na22_Ne22_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 12, 1, 11, 1, 6> j_Na22_Ne22_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Na22_Ne22_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 12> j_Na22_Ne22_temp;

    AMREX_GPU_MANAGED table_t j_Ne22_Na22_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 12, 1, 11, 1, 6> j_Ne22_Na22_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ne22_Na22_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 12> j_Ne22_Na22_temp;

    AMREX_GPU_MANAGED table_t j_Mg23_Na23_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 12, 1, 11, 1, 6> j_Mg23_Na23_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Mg23_Na23_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 12> j_Mg23_Na23_temp;

    AMREX_GPU_MANAGED table_t j_Na23_Mg23_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 12, 1, 11, 1, 6> j_Na23_Mg23_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Na23_Mg23_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 12> j_Na23_Mg23_temp;

    AMREX_GPU_MANAGED table_t j_Al25_Mg25_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 12, 1, 11, 1, 6> j_Al25_Mg25_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Al25_Mg25_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 12> j_Al25_Mg25_temp;

    AMREX_GPU_MANAGED table_t j_Mg25_Al25_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 12, 1, 11, 1, 6> j_Mg25_Al25_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Mg25_Al25_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 12> j_Mg25_Al25_temp;

    AMREX_GPU_MANAGED table_t j_Al26_Mg26_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 12, 1, 11, 1, 6> j_Al26_Mg26_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Al26_Mg26_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 12> j_Al26_Mg26_temp;

    AMREX_GPU_MANAGED table_t j_Mg26_Al26_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 12, 1, 11, 1, 6> j_Mg26_Al26_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Mg26_Al26_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 12> j_Mg26_Al26_temp;

    AMREX_GPU_MANAGED table_t j_Al27_Si27_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 12, 1, 11, 1, 6> j_Al27_Si27_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Al27_Si27_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 12> j_Al27_Si27_temp;

    AMREX_GPU_MANAGED table_t j_Si27_Al27_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 12, 1, 11, 1, 6> j_Si27_Al27_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Si27_Al27_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 12> j_Si27_Al27_temp;


}


void init_tabular()
{

    amrex::Print() << "reading in network electron-capture / beta-decay tables..." << std::endl;

    using namespace rate_tables;

    j_F17_O17_meta.ntemp = 12;
    j_F17_O17_meta.nrhoy = 11;
    j_F17_O17_meta.nvars = 6;
    j_F17_O17_meta.nheader = 5;

    init_tab_info(j_F17_O17_meta, "oda-17f-17o_electroncapture.dat", j_F17_O17_rhoy, j_F17_O17_temp, j_F17_O17_data);


    j_O17_F17_meta.ntemp = 12;
    j_O17_F17_meta.nrhoy = 11;
    j_O17_F17_meta.nvars = 6;
    j_O17_F17_meta.nheader = 5;

    init_tab_info(j_O17_F17_meta, "oda-17o-17f_betadecay.dat", j_O17_F17_rhoy, j_O17_F17_temp, j_O17_F17_data);


    j_F18_Ne18_meta.ntemp = 12;
    j_F18_Ne18_meta.nrhoy = 11;
    j_F18_Ne18_meta.nvars = 6;
    j_F18_Ne18_meta.nheader = 5;

    init_tab_info(j_F18_Ne18_meta, "oda-18f-18ne_betadecay.dat", j_F18_Ne18_rhoy, j_F18_Ne18_temp, j_F18_Ne18_data);


    j_F18_O18_meta.ntemp = 12;
    j_F18_O18_meta.nrhoy = 11;
    j_F18_O18_meta.nvars = 6;
    j_F18_O18_meta.nheader = 5;

    init_tab_info(j_F18_O18_meta, "oda-18f-18o_electroncapture.dat", j_F18_O18_rhoy, j_F18_O18_temp, j_F18_O18_data);


    j_Ne18_F18_meta.ntemp = 12;
    j_Ne18_F18_meta.nrhoy = 11;
    j_Ne18_F18_meta.nvars = 6;
    j_Ne18_F18_meta.nheader = 5;

    init_tab_info(j_Ne18_F18_meta, "oda-18ne-18f_electroncapture.dat", j_Ne18_F18_rhoy, j_Ne18_F18_temp, j_Ne18_F18_data);


    j_O18_F18_meta.ntemp = 12;
    j_O18_F18_meta.nrhoy = 11;
    j_O18_F18_meta.nvars = 6;
    j_O18_F18_meta.nheader = 5;

    init_tab_info(j_O18_F18_meta, "oda-18o-18f_betadecay.dat", j_O18_F18_rhoy, j_O18_F18_temp, j_O18_F18_data);


    j_F19_Ne19_meta.ntemp = 12;
    j_F19_Ne19_meta.nrhoy = 11;
    j_F19_Ne19_meta.nvars = 6;
    j_F19_Ne19_meta.nheader = 5;

    init_tab_info(j_F19_Ne19_meta, "oda-19f-19ne_betadecay.dat", j_F19_Ne19_rhoy, j_F19_Ne19_temp, j_F19_Ne19_data);


    j_Ne19_F19_meta.ntemp = 12;
    j_Ne19_F19_meta.nrhoy = 11;
    j_Ne19_F19_meta.nvars = 6;
    j_Ne19_F19_meta.nheader = 5;

    init_tab_info(j_Ne19_F19_meta, "oda-19ne-19f_electroncapture.dat", j_Ne19_F19_rhoy, j_Ne19_F19_temp, j_Ne19_F19_data);


    j_Na21_Ne21_meta.ntemp = 12;
    j_Na21_Ne21_meta.nrhoy = 11;
    j_Na21_Ne21_meta.nvars = 6;
    j_Na21_Ne21_meta.nheader = 5;

    init_tab_info(j_Na21_Ne21_meta, "oda-21na-21ne_electroncapture.dat", j_Na21_Ne21_rhoy, j_Na21_Ne21_temp, j_Na21_Ne21_data);


    j_Ne21_Na21_meta.ntemp = 12;
    j_Ne21_Na21_meta.nrhoy = 11;
    j_Ne21_Na21_meta.nvars = 6;
    j_Ne21_Na21_meta.nheader = 5;

    init_tab_info(j_Ne21_Na21_meta, "oda-21ne-21na_betadecay.dat", j_Ne21_Na21_rhoy, j_Ne21_Na21_temp, j_Ne21_Na21_data);


    j_Mg22_Na22_meta.ntemp = 12;
    j_Mg22_Na22_meta.nrhoy = 11;
    j_Mg22_Na22_meta.nvars = 6;
    j_Mg22_Na22_meta.nheader = 5;

    init_tab_info(j_Mg22_Na22_meta, "oda-22mg-22na_electroncapture.dat", j_Mg22_Na22_rhoy, j_Mg22_Na22_temp, j_Mg22_Na22_data);


    j_Na22_Mg22_meta.ntemp = 12;
    j_Na22_Mg22_meta.nrhoy = 11;
    j_Na22_Mg22_meta.nvars = 6;
    j_Na22_Mg22_meta.nheader = 5;

    init_tab_info(j_Na22_Mg22_meta, "oda-22na-22mg_betadecay.dat", j_Na22_Mg22_rhoy, j_Na22_Mg22_temp, j_Na22_Mg22_data);


    j_Na22_Ne22_meta.ntemp = 12;
    j_Na22_Ne22_meta.nrhoy = 11;
    j_Na22_Ne22_meta.nvars = 6;
    j_Na22_Ne22_meta.nheader = 5;

    init_tab_info(j_Na22_Ne22_meta, "oda-22na-22ne_electroncapture.dat", j_Na22_Ne22_rhoy, j_Na22_Ne22_temp, j_Na22_Ne22_data);


    j_Ne22_Na22_meta.ntemp = 12;
    j_Ne22_Na22_meta.nrhoy = 11;
    j_Ne22_Na22_meta.nvars = 6;
    j_Ne22_Na22_meta.nheader = 5;

    init_tab_info(j_Ne22_Na22_meta, "oda-22ne-22na_betadecay.dat", j_Ne22_Na22_rhoy, j_Ne22_Na22_temp, j_Ne22_Na22_data);


    j_Mg23_Na23_meta.ntemp = 12;
    j_Mg23_Na23_meta.nrhoy = 11;
    j_Mg23_Na23_meta.nvars = 6;
    j_Mg23_Na23_meta.nheader = 5;

    init_tab_info(j_Mg23_Na23_meta, "oda-23mg-23na_electroncapture.dat", j_Mg23_Na23_rhoy, j_Mg23_Na23_temp, j_Mg23_Na23_data);


    j_Na23_Mg23_meta.ntemp = 12;
    j_Na23_Mg23_meta.nrhoy = 11;
    j_Na23_Mg23_meta.nvars = 6;
    j_Na23_Mg23_meta.nheader = 5;

    init_tab_info(j_Na23_Mg23_meta, "oda-23na-23mg_betadecay.dat", j_Na23_Mg23_rhoy, j_Na23_Mg23_temp, j_Na23_Mg23_data);


    j_Al25_Mg25_meta.ntemp = 12;
    j_Al25_Mg25_meta.nrhoy = 11;
    j_Al25_Mg25_meta.nvars = 6;
    j_Al25_Mg25_meta.nheader = 5;

    init_tab_info(j_Al25_Mg25_meta, "oda-25al-25mg_electroncapture.dat", j_Al25_Mg25_rhoy, j_Al25_Mg25_temp, j_Al25_Mg25_data);


    j_Mg25_Al25_meta.ntemp = 12;
    j_Mg25_Al25_meta.nrhoy = 11;
    j_Mg25_Al25_meta.nvars = 6;
    j_Mg25_Al25_meta.nheader = 5;

    init_tab_info(j_Mg25_Al25_meta, "oda-25mg-25al_betadecay.dat", j_Mg25_Al25_rhoy, j_Mg25_Al25_temp, j_Mg25_Al25_data);


    j_Al26_Mg26_meta.ntemp = 12;
    j_Al26_Mg26_meta.nrhoy = 11;
    j_Al26_Mg26_meta.nvars = 6;
    j_Al26_Mg26_meta.nheader = 5;

    init_tab_info(j_Al26_Mg26_meta, "oda-26al-26mg_electroncapture.dat", j_Al26_Mg26_rhoy, j_Al26_Mg26_temp, j_Al26_Mg26_data);


    j_Mg26_Al26_meta.ntemp = 12;
    j_Mg26_Al26_meta.nrhoy = 11;
    j_Mg26_Al26_meta.nvars = 6;
    j_Mg26_Al26_meta.nheader = 5;

    init_tab_info(j_Mg26_Al26_meta, "oda-26mg-26al_betadecay.dat", j_Mg26_Al26_rhoy, j_Mg26_Al26_temp, j_Mg26_Al26_data);


    j_Al27_Si27_meta.ntemp = 12;
    j_Al27_Si27_meta.nrhoy = 11;
    j_Al27_Si27_meta.nvars = 6;
    j_Al27_Si27_meta.nheader = 5;

    init_tab_info(j_Al27_Si27_meta, "oda-27al-27si_betadecay.dat", j_Al27_Si27_rhoy, j_Al27_Si27_temp, j_Al27_Si27_data);


    j_Si27_Al27_meta.ntemp = 12;
    j_Si27_Al27_meta.nrhoy = 11;
    j_Si27_Al27_meta.nvars = 6;
    j_Si27_Al27_meta.nheader = 5;

    init_tab_info(j_Si27_Al27_meta, "oda-27si-27al_electroncapture.dat", j_Si27_Al27_rhoy, j_Si27_Al27_temp, j_Si27_Al27_data);



}
