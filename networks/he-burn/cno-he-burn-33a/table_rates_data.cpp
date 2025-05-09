#include <AMReX_Array.H>
#include <string>
#include <table_rates.H>
#include <AMReX_Print.H>

using namespace amrex;

namespace rate_tables
{

    AMREX_GPU_MANAGED table_t j_Mg22_Na22_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Mg22_Na22_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Mg22_Na22_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Mg22_Na22_temp;

    AMREX_GPU_MANAGED table_t j_Na22_Mg22_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 13, 1, 11, 1, 6> j_Na22_Mg22_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 11> j_Na22_Mg22_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 13> j_Na22_Mg22_temp;


}


void init_tabular()
{

    amrex::Print() << "reading in network electron-capture / beta-decay tables..." << std::endl;

    using namespace rate_tables;

    j_Mg22_Na22_meta.ntemp = 13;
    j_Mg22_Na22_meta.nrhoy = 11;
    j_Mg22_Na22_meta.nvars = 6;
    j_Mg22_Na22_meta.nheader = 5;

    init_tab_info(j_Mg22_Na22_meta, "ffn-22mg-22na_electroncapture.dat", j_Mg22_Na22_rhoy, j_Mg22_Na22_temp, j_Mg22_Na22_data);


    j_Na22_Mg22_meta.ntemp = 13;
    j_Na22_Mg22_meta.nrhoy = 11;
    j_Na22_Mg22_meta.nvars = 6;
    j_Na22_Mg22_meta.nheader = 5;

    init_tab_info(j_Na22_Mg22_meta, "ffn-22na-22mg_betadecay.dat", j_Na22_Mg22_rhoy, j_Na22_Mg22_temp, j_Na22_Mg22_data);



}
