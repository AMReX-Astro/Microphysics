#include <AMReX_Array.H>
#include <string>
#include <table_rates.H>
#include <AMReX_Print.H>

using namespace amrex;

namespace rate_tables
{

    AMREX_GPU_MANAGED table_t j_Ar39_K39_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Ar39_K39_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Ar39_K39_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Ar39_K39_temp;

    AMREX_GPU_MANAGED table_t j_K39_Ar39_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_K39_Ar39_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_K39_Ar39_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_K39_Ar39_temp;

    AMREX_GPU_MANAGED table_t j_P30_Si30_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_P30_Si30_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_P30_Si30_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_P30_Si30_temp;

    AMREX_GPU_MANAGED table_t j_P31_S31_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_P31_S31_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_P31_S31_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_P31_S31_temp;

    AMREX_GPU_MANAGED table_t j_S31_P31_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_S31_P31_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_S31_P31_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_S31_P31_temp;

    AMREX_GPU_MANAGED table_t j_Si30_P30_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Si30_P30_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Si30_P30_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Si30_P30_temp;

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


    j_P30_Si30_meta.ntemp = 13;
    j_P30_Si30_meta.nrhoy = 11;
    j_P30_Si30_meta.nvars = 6;
    j_P30_Si30_meta.nheader = 5;

    init_tab_info(j_P30_Si30_meta, "ffn-30p-30si_electroncapture.dat", j_P30_Si30_rhoy, j_P30_Si30_temp, j_P30_Si30_data);


    j_P31_S31_meta.ntemp = 13;
    j_P31_S31_meta.nrhoy = 11;
    j_P31_S31_meta.nvars = 6;
    j_P31_S31_meta.nheader = 5;

    init_tab_info(j_P31_S31_meta, "ffn-31p-31s_betadecay.dat", j_P31_S31_rhoy, j_P31_S31_temp, j_P31_S31_data);


    j_S31_P31_meta.ntemp = 13;
    j_S31_P31_meta.nrhoy = 11;
    j_S31_P31_meta.nvars = 6;
    j_S31_P31_meta.nheader = 5;

    init_tab_info(j_S31_P31_meta, "ffn-31s-31p_electroncapture.dat", j_S31_P31_rhoy, j_S31_P31_temp, j_S31_P31_data);


    j_Si30_P30_meta.ntemp = 13;
    j_Si30_P30_meta.nrhoy = 11;
    j_Si30_P30_meta.nvars = 6;
    j_Si30_P30_meta.nheader = 5;

    init_tab_info(j_Si30_P30_meta, "ffn-30si-30p_betadecay.dat", j_Si30_P30_rhoy, j_Si30_P30_temp, j_Si30_P30_data);


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
