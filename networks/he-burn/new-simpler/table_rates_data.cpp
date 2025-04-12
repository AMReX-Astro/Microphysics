#include <AMReX_Array.H>
#include <string>
#include <table_rates.H>
#include <AMReX_Print.H>

using namespace amrex;

namespace rate_tables
{

    AMREX_GPU_MANAGED table_t j_Mg23_Na23_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Mg23_Na23_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Mg23_Na23_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Mg23_Na23_temp;

    AMREX_GPU_MANAGED table_t j_Na22_Ne22_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Na22_Ne22_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Na22_Ne22_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Na22_Ne22_temp;

    AMREX_GPU_MANAGED table_t j_Na23_Mg23_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Na23_Mg23_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Na23_Mg23_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Na23_Mg23_temp;

    AMREX_GPU_MANAGED table_t j_Ne22_Na22_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ne22_Na22_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ne22_Na22_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ne22_Na22_temp;

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

    AMREX_GPU_MANAGED table_t j_n_p_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_n_p_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_n_p_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_n_p_temp;

    AMREX_GPU_MANAGED table_t j_Ni56_Co56_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ni56_Co56_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ni56_Co56_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ni56_Co56_temp;

    AMREX_GPU_MANAGED table_t j_p_n_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_p_n_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_p_n_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_p_n_temp;


}


void init_tabular()
{

    amrex::Print() << "reading in network electron-capture / beta-decay tables..." << std::endl;

    using namespace rate_tables;

    j_Mg23_Na23_meta.ntemp = 13;
    j_Mg23_Na23_meta.nrhoy = 11;
    j_Mg23_Na23_meta.nvars = 6;
    j_Mg23_Na23_meta.nheader = 5;

    init_tab_info(j_Mg23_Na23_meta, "ffn-23mg-23na_electroncapture.dat", j_Mg23_Na23_rhoy, j_Mg23_Na23_temp, j_Mg23_Na23_data);


    j_Na22_Ne22_meta.ntemp = 13;
    j_Na22_Ne22_meta.nrhoy = 11;
    j_Na22_Ne22_meta.nvars = 6;
    j_Na22_Ne22_meta.nheader = 5;

    init_tab_info(j_Na22_Ne22_meta, "ffn-22na-22ne_electroncapture.dat", j_Na22_Ne22_rhoy, j_Na22_Ne22_temp, j_Na22_Ne22_data);


    j_Na23_Mg23_meta.ntemp = 13;
    j_Na23_Mg23_meta.nrhoy = 11;
    j_Na23_Mg23_meta.nvars = 6;
    j_Na23_Mg23_meta.nheader = 5;

    init_tab_info(j_Na23_Mg23_meta, "ffn-23na-23mg_betadecay.dat", j_Na23_Mg23_rhoy, j_Na23_Mg23_temp, j_Na23_Mg23_data);


    j_Ne22_Na22_meta.ntemp = 13;
    j_Ne22_Na22_meta.nrhoy = 11;
    j_Ne22_Na22_meta.nvars = 6;
    j_Ne22_Na22_meta.nheader = 5;

    init_tab_info(j_Ne22_Na22_meta, "ffn-22ne-22na_betadecay.dat", j_Ne22_Na22_rhoy, j_Ne22_Na22_temp, j_Ne22_Na22_data);


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


    j_n_p_meta.ntemp = 13;
    j_n_p_meta.nrhoy = 11;
    j_n_p_meta.nvars = 6;
    j_n_p_meta.nheader = 5;

    init_tab_info(j_n_p_meta, "langanke-n-p_betadecay.dat", j_n_p_rhoy, j_n_p_temp, j_n_p_data);


    j_Ni56_Co56_meta.ntemp = 13;
    j_Ni56_Co56_meta.nrhoy = 11;
    j_Ni56_Co56_meta.nvars = 6;
    j_Ni56_Co56_meta.nheader = 5;

    init_tab_info(j_Ni56_Co56_meta, "langanke-56ni-56co_electroncapture.dat", j_Ni56_Co56_rhoy, j_Ni56_Co56_temp, j_Ni56_Co56_data);


    j_p_n_meta.ntemp = 13;
    j_p_n_meta.nrhoy = 11;
    j_p_n_meta.nvars = 6;
    j_p_n_meta.nheader = 5;

    init_tab_info(j_p_n_meta, "langanke-p-n_electroncapture.dat", j_p_n_rhoy, j_p_n_temp, j_p_n_data);



}
