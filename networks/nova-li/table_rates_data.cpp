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

    AMREX_GPU_MANAGED table_t j_F18_O18_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 12, 1, 11, 1, 6> j_F18_O18_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_F18_O18_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 12> j_F18_O18_temp;

    AMREX_GPU_MANAGED table_t j_Ne18_F18_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 12, 1, 11, 1, 6> j_Ne18_F18_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ne18_F18_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 12> j_Ne18_F18_temp;

    AMREX_GPU_MANAGED table_t j_Ne19_F19_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 12, 1, 11, 1, 6> j_Ne19_F19_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ne19_F19_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 12> j_Ne19_F19_temp;


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


    j_Ne19_F19_meta.ntemp = 12;
    j_Ne19_F19_meta.nrhoy = 11;
    j_Ne19_F19_meta.nvars = 6;
    j_Ne19_F19_meta.nheader = 5;

    init_tab_info(j_Ne19_F19_meta, "oda-19ne-19f_electroncapture.dat", j_Ne19_F19_rhoy, j_Ne19_F19_temp, j_Ne19_F19_data);



}
