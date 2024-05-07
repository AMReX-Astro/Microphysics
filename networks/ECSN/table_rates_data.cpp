#include <AMReX_Array.H>
#include <string>
#include <table_rates.H>
#include <AMReX_Print.H>

using namespace amrex;

namespace rate_tables
{

    AMREX_GPU_MANAGED table_t j_F20_O20_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 39, 1, 152, 1, 6> j_F20_O20_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 152> j_F20_O20_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 39> j_F20_O20_temp;

    AMREX_GPU_MANAGED table_t j_Ne20_F20_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 39, 1, 152, 1, 6> j_Ne20_F20_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 152> j_Ne20_F20_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 39> j_Ne20_F20_temp;

    AMREX_GPU_MANAGED table_t j_O20_F20_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 39, 1, 152, 1, 6> j_O20_F20_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 152> j_O20_F20_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 39> j_O20_F20_temp;

    AMREX_GPU_MANAGED table_t j_F20_Ne20_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 39, 1, 152, 1, 6> j_F20_Ne20_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 152> j_F20_Ne20_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 39> j_F20_Ne20_temp;


}


void init_tabular()
{

    amrex::Print() << "reading in network electron-capture / beta-decay tables..." << std::endl;

    using namespace rate_tables;

    j_F20_O20_meta.ntemp = 39;
    j_F20_O20_meta.nrhoy = 152;
    j_F20_O20_meta.nvars = 6;
    j_F20_O20_meta.nheader = 5;

    init_tab_info(j_F20_O20_meta, "20f-20o_electroncapture.dat", j_F20_O20_rhoy, j_F20_O20_temp, j_F20_O20_data);


    j_Ne20_F20_meta.ntemp = 39;
    j_Ne20_F20_meta.nrhoy = 152;
    j_Ne20_F20_meta.nvars = 6;
    j_Ne20_F20_meta.nheader = 7;

    init_tab_info(j_Ne20_F20_meta, "20ne-20f_electroncapture.dat", j_Ne20_F20_rhoy, j_Ne20_F20_temp, j_Ne20_F20_data);


    j_O20_F20_meta.ntemp = 39;
    j_O20_F20_meta.nrhoy = 152;
    j_O20_F20_meta.nvars = 6;
    j_O20_F20_meta.nheader = 6;

    init_tab_info(j_O20_F20_meta, "20o-20f_betadecay.dat", j_O20_F20_rhoy, j_O20_F20_temp, j_O20_F20_data);


    j_F20_Ne20_meta.ntemp = 39;
    j_F20_Ne20_meta.nrhoy = 152;
    j_F20_Ne20_meta.nvars = 6;
    j_F20_Ne20_meta.nheader = 7;

    init_tab_info(j_F20_Ne20_meta, "20f-20ne_betadecay.dat", j_F20_Ne20_rhoy, j_F20_Ne20_temp, j_F20_Ne20_data);



}
