#include <AMReX_Array.H>
#include <string>
#include <table_rates.H>
#include <AMReX_Print.H>

using namespace amrex;

namespace rate_tables
{

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
