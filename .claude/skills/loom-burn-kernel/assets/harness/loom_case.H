// Write a standalone case file (loom_case_format.H) from test_react: the
// burn inputs exactly as the Loom launcher would pass them (BurnParams,
// helmholtz tables, the initial state box) plus the stock Microphysics
// result of the same burn. The standalone reproducer (in the
// loom-burn-kernel-tools follow-up) reads it and needs nothing else from
// Microphysics or AMReX.
//
// Runtime parameter: unit_test.loom_dump_case=<file> (empty = off). Built
// whenever LOOM_KERNEL_DIR is set (USE_LOOM_VALIDATE=TRUE on the CPU, or
// USE_HIP=TRUE USE_LOOM_BURN=TRUE): the CPU build stores a stock CPU
// reference, the HIP build a stock HIP reference (the one the standalone
// hipcc kernel must match bit for bit). The initial state differs between
// the two as well: main.cpp builds the T/rho grid with std::pow/log10 in a
// ParallelFor, i.e. with glibc on the CPU and OCML on the GPU.
//
// Needs one box on one rank (n_cell = max_grid_size, no MPI), no ghost
// cells. The stock burn runs twice on copies of the state and the case is
// written only if both runs agree bit for bit; main.cpp's own burn then
// runs as usual on the untouched state.
#ifndef LOOM_CASE_H
#define LOOM_CASE_H

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_Gpu.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <AMReX_buildInfo.H>

#include <extern_parameters.H>
#include <burn_type.H>
#include <variables.H>
#include <react_zones.H>
#include <loom_params.H>
#include <loom_case_format.H>

// Stock burn of every zone: the do_react call of main.cpp's ParReduce body,
// with the per-zone success flag kept (zone_ok, 1 component).
inline void loom_case_stock_burn (amrex::MultiFab& state, amrex::iMultiFab& n_rhs_mf,
                                  amrex::iMultiFab& zone_ok, const plot_t& vars)
{
    for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox();
        amrex::Array4<amrex::Real> const sa = state.array(mfi);
        amrex::Array4<int> const na = n_rhs_mf.array(mfi);
        amrex::Array4<int> const oa = zone_ok.array(mfi);
        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            burn_t burn_state;
            const bool success = do_react(i, j, k, sa, burn_state, na, vars);
            oa(i, j, k) = success ? 1 : 0;
        });
    }
    amrex::Gpu::streamSynchronize();
}

// Host copy of the single local fab (no ghost cells: Array4 order, packed).
template <typename FAB>
std::vector<typename FAB::value_type> loom_case_to_host (const amrex::FabArray<FAB>& mf)
{
    std::vector<typename FAB::value_type> h;
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const FAB& fab = mf[mfi];
        h.resize(static_cast<std::size_t>(fab.size()));
        amrex::Gpu::dtoh_memcpy(h.data(), fab.dataPtr(), fab.nBytes());
    }
    return h;
}

inline void loom_case_write (std::FILE* fp, const void* p, std::size_t bytes, const std::string& path)
{
    if (bytes != 0 && std::fwrite(p, 1, bytes, fp) != bytes) {
        amrex::Abort("loom_dump_case: write to " + path + " failed");
    }
}

inline void loom_dump_case (const std::string& path, const amrex::MultiFab& state, const plot_t& vars)
{
    using namespace amrex;

    if (ParallelDescriptor::NProcs() > 1) {
        amrex::Abort("loom_dump_case: run on one MPI rank");
    }
    if (state.size() != 1 || state.nGrow() != 0) {
        amrex::Abort("loom_dump_case: needs exactly one box without ghost cells (max_grid_size = n_cell)");
    }

    const BoxArray& ba = state.boxArray();
    const DistributionMapping& dm = state.DistributionMap();
    const int ncomp = state.nComp();
    const Box bx = ba[0];
    const long ncells = static_cast<long>(bx.numPts());

    // the kernel's parameters and tables, exactly as the launcher fills them
    mpl::BurnParams hp{};
    fill_burn_params(&hp);        // aborts on unsupported settings
    mpl::EosTables tab{};
    fill_eos_tables(&tab);
    const mpl::Real* tables[LOOM_CASE_NTABLES] = {
        tab.f, tab.dpdf, tab.ef, tab.xf, tab.d, tab.t, tab.dt_sav, tab.dt2_sav,
        tab.dti_sav, tab.dt2i_sav, tab.dd_sav, tab.dd2_sav, tab.ddi_sav, tab.dd2i_sav};
    const std::uint64_t table_len[LOOM_CASE_NTABLES] = {
        sizeof(helmholtz::f) / sizeof(Real), sizeof(helmholtz::dpdf) / sizeof(Real),
        sizeof(helmholtz::ef) / sizeof(Real), sizeof(helmholtz::xf) / sizeof(Real),
        sizeof(helmholtz::d) / sizeof(Real), sizeof(helmholtz::t) / sizeof(Real),
        sizeof(helmholtz::dt_sav) / sizeof(Real), sizeof(helmholtz::dt2_sav) / sizeof(Real),
        sizeof(helmholtz::dti_sav) / sizeof(Real), sizeof(helmholtz::dt2i_sav) / sizeof(Real),
        sizeof(helmholtz::dd_sav) / sizeof(Real), sizeof(helmholtz::dd2_sav) / sizeof(Real),
        sizeof(helmholtz::ddi_sav) / sizeof(Real), sizeof(helmholtz::dd2i_sav) / sizeof(Real)};

    // two stock burns of copies of the initial state; they must agree
    std::vector<Real> ref_state;
    std::vector<int> ref_n_rhs;
    std::vector<int> ref_ok;
    for (int run = 0; run < 2; ++run) {
        MultiFab s(ba, dm, ncomp, 0);
        MultiFab::Copy(s, state, 0, 0, ncomp, 0);
        iMultiFab n(ba, dm, 2, 0);
        n.setVal(-1);
        iMultiFab ok(ba, dm, 1, 0);
        ok.setVal(-1);
        loom_case_stock_burn(s, n, ok, vars);
        std::vector<Real> hs = loom_case_to_host(s);
        std::vector<int> hn = loom_case_to_host(n);
        std::vector<int> ho = loom_case_to_host(ok);
        if (run == 0) {
            ref_state = std::move(hs);
            ref_n_rhs = std::move(hn);
            ref_ok = std::move(ho);
        } else if (std::memcmp(hs.data(), ref_state.data(), hs.size() * sizeof(Real)) != 0 ||
                   hn != ref_n_rhs || ho != ref_ok) {
            amrex::Abort("loom_dump_case: two stock burns of the same state differ; not writing " + path);
        }
    }
    const std::vector<Real> state_in = loom_case_to_host(state);
    int ref_failures = 0;
    for (int v : ref_ok) {
        if (v != 1) { ++ref_failures; }
    }

    LoomCaseHeader h{};
    std::memcpy(h.magic, LOOM_CASE_MAGIC, sizeof h.magic);
    h.version = LOOM_CASE_VERSION;
    h.header_bytes = sizeof(LoomCaseHeader);
    h.num_spec = NumSpec;
    h.ncomp = ncomp;
    h.nx = bx.length(0);
    h.ny = bx.length(1);
    h.nz = bx.length(2);
    h.lo_i = bx.smallEnd(0);
    h.lo_j = bx.smallEnd(1);
    h.lo_k = bx.smallEnd(2);
    h.irho = vars.irho;
    h.itemp = vars.itemp;
    h.ispec_old = vars.ispec_old;
    h.ispec = vars.ispec;
    h.irodot = vars.irodot;
    h.irho_hnuc = vars.irho_hnuc;
#ifdef AMREX_USE_GPU
    h.ref_kind = LOOM_CASE_REF_STOCK_HIP;
#else
    h.ref_kind = LOOM_CASE_REF_STOCK_CPU;
#endif
    h.ref_failures = ref_failures;
    h.dt = unit_test_rp::tmax;
    h.burn_params_bytes = sizeof(mpl::BurnParams);
    for (int t = 0; t < LOOM_CASE_NTABLES; ++t) { h.table_len[t] = table_len[t]; }

    std::string desc = "stock Microphysics test_react";
#ifdef AMREX_USE_GPU
    desc += " (HIP)";
#else
    desc += " (CPU)";
#endif
    for (int i = 1; i <= buildInfoGetNumModules(); ++i) {
        desc += std::string(" ") + buildInfoGetModuleName(i) + "=" + buildInfoGetModuleVal(i);
    }
    desc += std::string(" kernel_dir=") + LOOM_KERNEL_STR(LOOM_KERNEL_DIR_PATH);
    for (int i = 1; i <= 2; ++i) {   // empty when the tree has no git metadata
        const std::string g = buildInfoGetGitHash(i);
        if (!g.empty()) { desc += " git" + std::to_string(i) + "=" + g; }
    }
    desc += std::string(" comp=") + buildInfoGetComp() + " " + buildInfoGetCompVersion();
#ifdef AMREX_USE_HIP
    {
        hipDeviceProp_t prop{};
        int dev = 0;
        if (hipGetDevice(&dev) == hipSuccess && hipGetDeviceProperties(&prop, dev) == hipSuccess) {
            desc += std::string(" device=") + prop.gcnArchName;
        }
    }
#endif
    desc += std::string(" cxxflags=") + buildInfoGetCXXFlags();
    std::strncpy(h.description, desc.c_str(), sizeof h.description - 1);

    std::FILE* fp = std::fopen(path.c_str(), "wb");
    if (fp == nullptr) {
        amrex::Abort("loom_dump_case: cannot open " + path);
    }
    loom_case_write(fp, &h, sizeof h, path);
    loom_case_write(fp, &hp, sizeof hp, path);
    for (int t = 0; t < LOOM_CASE_NTABLES; ++t) {
        loom_case_write(fp, tables[t], table_len[t] * sizeof(Real), path);
    }
    loom_case_write(fp, state_in.data(), state_in.size() * sizeof(Real), path);
    loom_case_write(fp, ref_state.data(), ref_state.size() * sizeof(Real), path);
    loom_case_write(fp, ref_n_rhs.data(), ref_n_rhs.size() * sizeof(int), path);
    loom_case_write(fp, ref_ok.data(), ref_ok.size() * sizeof(int), path);
    if (std::fclose(fp) != 0) {
        amrex::Abort("loom_dump_case: closing " + path + " failed");
    }

    amrex::Print() << "LOOM_DUMP_CASE file=" << path << " zones=" << ncells << " ncomp=" << ncomp
                   << " ref_failures=" << ref_failures << " burn_params_bytes=" << sizeof(mpl::BurnParams)
                   << " (" << desc.substr(0, 120) << " ...)\n";
}

#endif
