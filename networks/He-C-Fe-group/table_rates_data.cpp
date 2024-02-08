#include <AMReX_Array.H>
#include <string>
#include <table_rates.H>
#include <AMReX_Print.H>

using namespace amrex;

namespace rate_tables
{

    AMREX_GPU_MANAGED table_t j_Co55_Fe55_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Co55_Fe55_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Co55_Fe55_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Co55_Fe55_temp;

    AMREX_GPU_MANAGED table_t j_Co56_Fe56_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Co56_Fe56_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Co56_Fe56_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Co56_Fe56_temp;

    AMREX_GPU_MANAGED table_t j_Co56_Ni56_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Co56_Ni56_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Co56_Ni56_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Co56_Ni56_temp;

    AMREX_GPU_MANAGED table_t j_Co57_Ni57_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Co57_Ni57_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Co57_Ni57_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Co57_Ni57_temp;

    AMREX_GPU_MANAGED table_t j_Fe55_Co55_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Fe55_Co55_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Fe55_Co55_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Fe55_Co55_temp;

    AMREX_GPU_MANAGED table_t j_Fe55_Mn55_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Fe55_Mn55_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Fe55_Mn55_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Fe55_Mn55_temp;

    AMREX_GPU_MANAGED table_t j_Fe56_Co56_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Fe56_Co56_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Fe56_Co56_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Fe56_Co56_temp;

    AMREX_GPU_MANAGED table_t j_Mn55_Fe55_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Mn55_Fe55_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Mn55_Fe55_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Mn55_Fe55_temp;

    AMREX_GPU_MANAGED table_t j_n_p_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_n_p_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_n_p_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_n_p_temp;

    AMREX_GPU_MANAGED table_t j_Ni56_Co56_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ni56_Co56_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ni56_Co56_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ni56_Co56_temp;

    AMREX_GPU_MANAGED table_t j_Ni57_Co57_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ni57_Co57_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ni57_Co57_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ni57_Co57_temp;

    AMREX_GPU_MANAGED table_t j_p_n_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_p_n_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_p_n_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_p_n_temp;


}


void init_tabular()
{

    amrex::Print() << "reading in network electron-capture / beta-decay tables..." << std::endl;

    using namespace rate_tables;

    j_Co55_Fe55_meta.ntemp = 13;
    j_Co55_Fe55_meta.nrhoy = 11;
    j_Co55_Fe55_meta.nvars = 6;
    j_Co55_Fe55_meta.nheader = 5;

    init_tab_info(j_Co55_Fe55_meta, "55co-55fe_electroncapture.dat", j_Co55_Fe55_rhoy, j_Co55_Fe55_temp, j_Co55_Fe55_data);


    j_Co56_Fe56_meta.ntemp = 13;
    j_Co56_Fe56_meta.nrhoy = 11;
    j_Co56_Fe56_meta.nvars = 6;
    j_Co56_Fe56_meta.nheader = 5;

    init_tab_info(j_Co56_Fe56_meta, "56co-56fe_electroncapture.dat", j_Co56_Fe56_rhoy, j_Co56_Fe56_temp, j_Co56_Fe56_data);


    j_Co56_Ni56_meta.ntemp = 13;
    j_Co56_Ni56_meta.nrhoy = 11;
    j_Co56_Ni56_meta.nvars = 6;
    j_Co56_Ni56_meta.nheader = 5;

    init_tab_info(j_Co56_Ni56_meta, "56co-56ni_betadecay.dat", j_Co56_Ni56_rhoy, j_Co56_Ni56_temp, j_Co56_Ni56_data);


    j_Co57_Ni57_meta.ntemp = 13;
    j_Co57_Ni57_meta.nrhoy = 11;
    j_Co57_Ni57_meta.nvars = 6;
    j_Co57_Ni57_meta.nheader = 5;

    init_tab_info(j_Co57_Ni57_meta, "57co-57ni_betadecay.dat", j_Co57_Ni57_rhoy, j_Co57_Ni57_temp, j_Co57_Ni57_data);


    j_Fe55_Co55_meta.ntemp = 13;
    j_Fe55_Co55_meta.nrhoy = 11;
    j_Fe55_Co55_meta.nvars = 6;
    j_Fe55_Co55_meta.nheader = 5;

    init_tab_info(j_Fe55_Co55_meta, "55fe-55co_betadecay.dat", j_Fe55_Co55_rhoy, j_Fe55_Co55_temp, j_Fe55_Co55_data);


    j_Fe55_Mn55_meta.ntemp = 13;
    j_Fe55_Mn55_meta.nrhoy = 11;
    j_Fe55_Mn55_meta.nvars = 6;
    j_Fe55_Mn55_meta.nheader = 5;

    init_tab_info(j_Fe55_Mn55_meta, "55fe-55mn_electroncapture.dat", j_Fe55_Mn55_rhoy, j_Fe55_Mn55_temp, j_Fe55_Mn55_data);


    j_Fe56_Co56_meta.ntemp = 13;
    j_Fe56_Co56_meta.nrhoy = 11;
    j_Fe56_Co56_meta.nvars = 6;
    j_Fe56_Co56_meta.nheader = 5;

    init_tab_info(j_Fe56_Co56_meta, "56fe-56co_betadecay.dat", j_Fe56_Co56_rhoy, j_Fe56_Co56_temp, j_Fe56_Co56_data);


    j_Mn55_Fe55_meta.ntemp = 13;
    j_Mn55_Fe55_meta.nrhoy = 11;
    j_Mn55_Fe55_meta.nvars = 6;
    j_Mn55_Fe55_meta.nheader = 5;

    init_tab_info(j_Mn55_Fe55_meta, "55mn-55fe_betadecay.dat", j_Mn55_Fe55_rhoy, j_Mn55_Fe55_temp, j_Mn55_Fe55_data);


    j_n_p_meta.ntemp = 13;
    j_n_p_meta.nrhoy = 11;
    j_n_p_meta.nvars = 6;
    j_n_p_meta.nheader = 5;

    init_tab_info(j_n_p_meta, "n-p_betadecay.dat", j_n_p_rhoy, j_n_p_temp, j_n_p_data);


    j_Ni56_Co56_meta.ntemp = 13;
    j_Ni56_Co56_meta.nrhoy = 11;
    j_Ni56_Co56_meta.nvars = 6;
    j_Ni56_Co56_meta.nheader = 5;

    init_tab_info(j_Ni56_Co56_meta, "56ni-56co_electroncapture.dat", j_Ni56_Co56_rhoy, j_Ni56_Co56_temp, j_Ni56_Co56_data);


    j_Ni57_Co57_meta.ntemp = 13;
    j_Ni57_Co57_meta.nrhoy = 11;
    j_Ni57_Co57_meta.nvars = 6;
    j_Ni57_Co57_meta.nheader = 5;

    init_tab_info(j_Ni57_Co57_meta, "57ni-57co_electroncapture.dat", j_Ni57_Co57_rhoy, j_Ni57_Co57_temp, j_Ni57_Co57_data);


    j_p_n_meta.ntemp = 13;
    j_p_n_meta.nrhoy = 11;
    j_p_n_meta.nvars = 6;
    j_p_n_meta.nheader = 5;

    init_tab_info(j_p_n_meta, "p-n_electroncapture.dat", j_p_n_rhoy, j_p_n_temp, j_p_n_data);



}
