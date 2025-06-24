#include <AMReX_Array.H>
#include <string>
#include <table_rates.H>
#include <AMReX_Print.H>

using namespace amrex;

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

    AMREX_GPU_MANAGED table_t j_Na22_Ne22_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 12, 1, 11, 1, 6> j_Na22_Ne22_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Na22_Ne22_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 12> j_Na22_Ne22_temp;

    AMREX_GPU_MANAGED table_t j_Ne22_Na22_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 12, 1, 11, 1, 6> j_Ne22_Na22_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ne22_Na22_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 12> j_Ne22_Na22_temp;


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



}
