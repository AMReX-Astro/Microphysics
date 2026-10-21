// HIP host launcher for the Loom burn kernel (loom/<network>_<integrator>/burn_kernel.cpp,
// root mpl::burn_zones). Compiled only with USE_HIP=TRUE USE_LOOM_BURN=TRUE
// (-DMICROPHYSICS_LOOM); see loom/<network>_<integrator>/README.md and
// references/test-react-integration.md in the loom-burn-kernel skill.
//
// Runtime parameters (namespace unit_test):
//   loom_hsaco    code object (.hsaco) from loom-compile; empty = stock path
//   loom_kernel   kernel symbol in the code object (default "mpl.burn_zones")
//   loom_ncells   zones per box the code object was compiled for
//                 (workgroup_count.x = ceil(loom_ncells / 128)); every box
//                 must have exactly this many zones (max_grid_size = n_cell)
//   loom_compare  1: run stock, Loom, stock, Loom, ... (max(2, loom_repeats)
//                 pairs) from the same initial state and compare bit for bit;
//                 prints LOOM_COMPARE and aborts unless everything is identical
//   loom_repeats  number of Loom burns (timing); pairs in compare mode
//
// Device data:
//   - BurnParams: filled on the host from the stock globals (loom_params.H),
//     copied to an amrex::The_Arena() allocation.
//   - helmholtz tables: in a HIP build the stock arrays are AMREX_GPU_MANAGED
//     (__managed__), so their addresses are valid on the device and are
//     passed directly. Without managed memory, copy each array (same layout)
//     to device memory and pass those pointers instead.
//   - state / n_rhs: the MultiFab / iMultiFab box data (Array4 pointer and
//     strides); n_rhs has 2 components (n_rhs, n_step).
//   - zone_success: one int per zone, reduced on the host into num_failed.
//
// zone_state.out: the zone with the most RHS evaluations is found on the host
// from n_rhs (first such zone in AMReX order, i fastest); its burn_t is
// recomputed by a HOST re-burn of that zone with the stock do_react (the
// kernel does not return burn_t), so its contents are the stock CPU result.
#ifndef LOOM_LAUNCH_H
#define LOOM_LAUNCH_H

#ifdef MICROPHYSICS_LOOM

#if !defined(AMREX_USE_HIP)
#error "MICROPHYSICS_LOOM (USE_LOOM_BURN=TRUE) requires USE_HIP=TRUE"
#endif

#include <hip/hip_runtime.h>

#include <cmath>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_Gpu.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <AMReX_ValLocPair.H>

#include <extern_parameters.H>
#include <burn_type.H>
#include <variables.H>
#include <react_zones.H>
#include <loom_params.H>

// The kernel root's host form (an inline function with one extra trailing
// int, the zone index): used only to check at compile time that the
// argument array below matches the kernel's declared parameter types.
// Reached by path from this directory (LOOM_KERNEL_PATH, defined in
// <loom_eos_params.H> via <loom_params.H> from LOOM_KERNEL_DIR on the make
// command line).
#include LOOM_KERNEL_PATH(burn_kernel.cpp)
#include <type_traits>

#define LOOM_HIP_CHECK(call)                                                        \
    do {                                                                            \
        hipError_t loom_err_ = (call);                                              \
        if (loom_err_ != hipSuccess) {                                              \
            amrex::Abort(std::string("Loom launcher: ") + #call + " failed: " +     \
                         hipGetErrorString(loom_err_));                             \
        }                                                                           \
    } while (0)

// One Loom burn of every local box.
struct LoomRunStats {
    int num_failed{0};
    float ms{0.0f};        // kernel time (hipEvents), summed over boxes
};

class LoomBurner
{
public:
    LoomBurner (const std::string& hsaco, const std::string& kernel, int ncells)
        : m_ncells(ncells)
    {
        if (ncells <= 0) {
            amrex::Abort("Loom launcher: set unit_test.loom_ncells to the zone count the code object was compiled for");
        }

        mpl::BurnParams hp{};
        fill_burn_params(&hp);    // aborts on unsupported settings (network.use_tables)
        m_params_d = static_cast<mpl::BurnParams*>(amrex::The_Arena()->alloc(sizeof(mpl::BurnParams)));
        amrex::Gpu::htod_memcpy(m_params_d, &hp, sizeof(mpl::BurnParams));

        // The table addresses are passed as is, which needs managed memory
        // on this device (the stock HIP burn reads the same arrays).
        int dev = 0, managed = 0;
        LOOM_HIP_CHECK(hipGetDevice(&dev));
        LOOM_HIP_CHECK(hipDeviceGetAttribute(&managed, hipDeviceAttributeManagedMemory, dev));
        if (managed == 0) {
            amrex::Abort("Loom launcher: device has no managed memory; the AMREX_GPU_MANAGED helmholtz tables would need device copies");
        }
        fill_eos_tables(&m_tab);  // managed arrays: valid device addresses

        m_success_d = static_cast<int*>(amrex::The_Arena()->alloc(sizeof(int) * static_cast<std::size_t>(ncells)));

        LOOM_HIP_CHECK(hipModuleLoad(&m_mod, hsaco.c_str()));
        LOOM_HIP_CHECK(hipModuleGetFunction(&m_fn, m_mod, kernel.c_str()));
        LOOM_HIP_CHECK(hipEventCreate(&m_start));
        LOOM_HIP_CHECK(hipEventCreate(&m_stop));

        amrex::Print() << "Loom launcher: loaded " << hsaco << " kernel " << kernel
                       << " ncells=" << ncells << " blocks=" << blocks() << " threads=128\n";
    }

    ~LoomBurner ()
    {
        (void) hipEventDestroy(m_start);
        (void) hipEventDestroy(m_stop);
        (void) hipModuleUnload(m_mod);
        amrex::The_Arena()->free(m_success_d);
        amrex::The_Arena()->free(m_params_d);
    }

    LoomBurner (const LoomBurner&) = delete;
    LoomBurner& operator= (const LoomBurner&) = delete;

    [[nodiscard]] unsigned int blocks () const { return static_cast<unsigned int>((m_ncells + 127) / 128); }

    LoomRunStats burn (amrex::MultiFab& state, amrex::iMultiFab& n_rhs_mf, const plot_t& vars, amrex::Real dt_in)
    {
        AMREX_ALWAYS_ASSERT(state.nGrow() == 0 && n_rhs_mf.nGrow() == 0 && n_rhs_mf.nComp() == 2);

        LoomRunStats rs;
        const hipStream_t stream = amrex::Gpu::gpuStream();
        std::vector<int> success_h(static_cast<std::size_t>(m_ncells));

        for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
            const amrex::Box bx = mfi.validbox();
            if (bx.numPts() != static_cast<amrex::Long>(m_ncells)) {
                amrex::Abort("Loom launcher: box " + std::to_string(bx.numPts()) +
                             " zones != unit_test.loom_ncells = " + std::to_string(m_ncells) +
                             " (run with max_grid_size = n_cell, or compile the code object for this box size)");
            }
            amrex::Array4<amrex::Real> const sa = state.array(mfi);
            amrex::Array4<int> const na = n_rhs_mf.array(mfi);
            AMREX_ALWAYS_ASSERT(sa.begin.vect[0] == bx.smallEnd(0) && sa.begin.vect[1] == bx.smallEnd(1) &&
                                sa.begin.vect[2] == bx.smallEnd(2));
            AMREX_ALWAYS_ASSERT(na.begin.vect[0] == bx.smallEnd(0) && na.begin.vect[1] == bx.smallEnd(1) &&
                                na.begin.vect[2] == bx.smallEnd(2));

            // Each argument in a variable of the exact type burn_zones declares.
            const mpl::BurnParams* a_params = m_params_d;
            const mpl::Real* a_f = m_tab.f;
            const mpl::Real* a_dpdf = m_tab.dpdf;
            const mpl::Real* a_ef = m_tab.ef;
            const mpl::Real* a_xf = m_tab.xf;
            const mpl::Real* a_d = m_tab.d;
            const mpl::Real* a_t = m_tab.t;
            const mpl::Real* a_dt_sav = m_tab.dt_sav;
            const mpl::Real* a_dt2_sav = m_tab.dt2_sav;
            const mpl::Real* a_dti_sav = m_tab.dti_sav;
            const mpl::Real* a_dt2i_sav = m_tab.dt2i_sav;
            const mpl::Real* a_dd_sav = m_tab.dd_sav;
            const mpl::Real* a_dd2_sav = m_tab.dd2_sav;
            const mpl::Real* a_ddi_sav = m_tab.ddi_sav;
            const mpl::Real* a_dd2i_sav = m_tab.dd2i_sav;
            mpl::Real* a_state = sa.p;
            long a_jstride = static_cast<long>(sa.get_stride<1>());
            long a_kstride = static_cast<long>(sa.get_stride<2>());
            long a_nstride = static_cast<long>(sa.get_stride<3>());
            int a_lo_i = bx.smallEnd(0);
            int a_lo_j = bx.smallEnd(1);
            int a_lo_k = bx.smallEnd(2);
            int a_nx = bx.length(0);
            int a_ny = bx.length(1);
            int a_nz = bx.length(2);
            int a_irho = vars.irho;
            int a_itemp = vars.itemp;
            int a_ispec_old = vars.ispec_old;
            int a_ispec = vars.ispec;
            int a_irodot = vars.irodot;
            int a_irho_hnuc = vars.irho_hnuc;
            int* a_n_rhs = na.p;
            long a_n_rhs_jstride = static_cast<long>(na.get_stride<1>());
            long a_n_rhs_kstride = static_cast<long>(na.get_stride<2>());
            long a_n_rhs_nstride = static_cast<long>(na.get_stride<3>());
            mpl::Real a_dt = dt_in;
            int* a_zone_success = m_success_d;

            // the argument types must be exactly the kernel's parameter types
            using KernelSig = void (*)(decltype(a_params),
                decltype(a_f), decltype(a_dpdf), decltype(a_ef), decltype(a_xf), decltype(a_d), decltype(a_t),
                decltype(a_dt_sav), decltype(a_dt2_sav), decltype(a_dti_sav), decltype(a_dt2i_sav),
                decltype(a_dd_sav), decltype(a_dd2_sav), decltype(a_ddi_sav), decltype(a_dd2i_sav),
                decltype(a_state), decltype(a_jstride), decltype(a_kstride), decltype(a_nstride),
                decltype(a_lo_i), decltype(a_lo_j), decltype(a_lo_k), decltype(a_nx), decltype(a_ny), decltype(a_nz),
                decltype(a_irho), decltype(a_itemp), decltype(a_ispec_old), decltype(a_ispec), decltype(a_irodot),
                decltype(a_irho_hnuc),
                decltype(a_n_rhs), decltype(a_n_rhs_jstride), decltype(a_n_rhs_kstride), decltype(a_n_rhs_nstride),
                decltype(a_dt), decltype(a_zone_success),
                int /* host-only zone index */);
            static_assert(std::is_same_v<decltype(&mpl::burn_zones), KernelSig>,
                          "Loom launcher arguments do not match mpl::burn_zones");

            void* args[] = {
                &a_params,
                &a_f, &a_dpdf, &a_ef, &a_xf, &a_d, &a_t, &a_dt_sav, &a_dt2_sav,
                &a_dti_sav, &a_dt2i_sav, &a_dd_sav, &a_dd2_sav, &a_ddi_sav, &a_dd2i_sav,
                &a_state, &a_jstride, &a_kstride, &a_nstride,
                &a_lo_i, &a_lo_j, &a_lo_k, &a_nx, &a_ny, &a_nz,
                &a_irho, &a_itemp, &a_ispec_old, &a_ispec, &a_irodot, &a_irho_hnuc,
                &a_n_rhs, &a_n_rhs_jstride, &a_n_rhs_kstride, &a_n_rhs_nstride,
                &a_dt, &a_zone_success
            };

            LOOM_HIP_CHECK(hipEventRecord(m_start, stream));
            LOOM_HIP_CHECK(hipModuleLaunchKernel(m_fn, blocks(), 1, 1, 128, 1, 1, 0, stream, args, nullptr));
            LOOM_HIP_CHECK(hipGetLastError());
            LOOM_HIP_CHECK(hipEventRecord(m_stop, stream));
            LOOM_HIP_CHECK(hipEventSynchronize(m_stop));
            float ms = 0.0f;
            LOOM_HIP_CHECK(hipEventElapsedTime(&ms, m_start, m_stop));
            rs.ms += ms;

            amrex::Gpu::dtoh_memcpy(success_h.data(), m_success_d, sizeof(int) * success_h.size());
            for (int s : success_h) {
                if (s == 0) { ++rs.num_failed; }
            }
        }
        amrex::Gpu::streamSynchronize();
        return rs;
    }

private:
    int m_ncells{0};
    hipModule_t m_mod{};
    hipFunction_t m_fn{};
    hipEvent_t m_start{};
    hipEvent_t m_stop{};
    mpl::BurnParams* m_params_d{nullptr};
    int* m_success_d{nullptr};
    mpl::EosTables m_tab{};
};

// The stock burn of main.cpp (same ParReduce body), for loom_compare.
inline amrex::ValLocPair<int, burn_t>
loom_stock_burn (amrex::MultiFab& state, amrex::iMultiFab& integrator_n_rhs, const plot_t& vars,
                 int& num_failed, float& ms)
{
    using namespace amrex;

    num_failed = 0;
    AsyncArray<int> aa_num_failed(&num_failed, 1);
    int* num_failed_d = aa_num_failed.data();

    hipEvent_t start{};
    hipEvent_t stop{};
    LOOM_HIP_CHECK(hipEventCreate(&start));
    LOOM_HIP_CHECK(hipEventCreate(&stop));
    LOOM_HIP_CHECK(hipEventRecord(start, Gpu::gpuStream()));

    ValLocPair<int, burn_t> r;
    {
        auto const& ma = state.arrays();
        auto const& ia = integrator_n_rhs.arrays();

        r = amrex::ParReduce(TypeList<ReduceOpMax>{}, TypeList<ValLocPair<int, burn_t>>{}, state,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) -> GpuTuple<ValLocPair<int, burn_t>>
        {
            Array4<Real> const& s = ma[box_no];
            auto n_rhs = ia[box_no];

            burn_t burn_state;
            bool success = do_react(i, j, k, s, burn_state, n_rhs, vars);

            if (!success) {
                Gpu::Atomic::Add(num_failed_d, 1);
            }

            return {ValLocPair<int, burn_t>{n_rhs(i,j,k,0), burn_state}};
        });
    }

    LOOM_HIP_CHECK(hipEventRecord(stop, Gpu::gpuStream()));
    LOOM_HIP_CHECK(hipEventSynchronize(stop));
    LOOM_HIP_CHECK(hipEventElapsedTime(&ms, start, stop));
    (void) hipEventDestroy(start);
    (void) hipEventDestroy(stop);

    aa_num_failed.copyToHost(&num_failed, 1);
    Gpu::synchronize();
    return r;
}

// Host copy of every local fab of a FabArray (box order of the MFIter).
template <typename FAB>
std::vector<std::vector<typename FAB::value_type>> loom_to_host (const amrex::FabArray<FAB>& mf)
{
    std::vector<std::vector<typename FAB::value_type>> h;
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const FAB& fab = mf[mfi];
        h.emplace_back(static_cast<std::size_t>(fab.size()));
        amrex::Gpu::dtoh_memcpy(h.back().data(), fab.dataPtr(), fab.nBytes());
    }
    return h;
}

// Zone with the most RHS evaluations (first in AMReX order), and the burn_t
// of a host stock re-burn of that zone (zone_state.out).
inline amrex::ValLocPair<int, burn_t>
loom_max_nrhs_zone (const amrex::MultiFab& state, const amrex::iMultiFab& integrator_n_rhs, const plot_t& vars)
{
    using namespace amrex;
    const auto nh = loom_to_host(integrator_n_rhs);
    int best = std::numeric_limits<int>::lowest();
    IntVect best_iv(0);
    int best_box = -1;
    int ib = 0;
    for (MFIter mfi(integrator_n_rhs); mfi.isValid(); ++mfi, ++ib) {
        const Box bx = mfi.validbox();
        const Dim3 lo = amrex::lbound(bx);
        const Dim3 hi{bx.bigEnd(0) + 1, bx.bigEnd(1) + 1, bx.bigEnd(2) + 1};
        Array4<int const> const a(nh[ib].data(), lo, hi, 2);
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
            if (a(i, j, k, 0) > best) { best = a(i, j, k, 0); best_iv = IntVect(i, j, k); best_box = ib; }
        });
    }

    ValLocPair<int, burn_t> r{best, burn_t{}};
    if (best_box < 0) { return r; }

    // one-zone host copy of the (unchanged) inputs, then the stock do_react
    const int ncomp = state.nComp();
    std::vector<Real> zone(static_cast<std::size_t>(ncomp));
    std::vector<int> zone_n(2);
    ib = 0;
    for (MFIter mfi(state); mfi.isValid(); ++mfi, ++ib) {
        if (ib != best_box) { continue; }
        Array4<Real const> const sa = state.const_array(mfi);
        for (int n = 0; n < ncomp; ++n) {
            Gpu::dtoh_memcpy(&zone[n], sa.ptr(best_iv[0], best_iv[1], best_iv[2], n), sizeof(Real));
        }
    }
    const Dim3 lo{best_iv[0], best_iv[1], best_iv[2]};
    const Dim3 hi{best_iv[0] + 1, best_iv[1] + 1, best_iv[2] + 1};
    Array4<Real> const za(zone.data(), lo, hi, ncomp);
    Array4<int> const zn(zone_n.data(), lo, hi, 2);
    burn_t burn_state;
    do_react(best_iv[0], best_iv[1], best_iv[2], za, burn_state, zn, vars);
    r.index = burn_state;
    return r;
}

// Per-component and per-zone comparison of two burns (stock vs Loom).
struct LoomCompStat {
    std::string name;
    double max_abs{0.0};
    double max_rel{0.0};
    int i{0}, j{0}, k{0};
    double rho{0.0}, T{0.0};
    long differing{0};
};

struct LoomCompare {
    long zones{0};
    long identical_zones{0};
    long nrhs_diff_zones{0};
    long nstep_diff_zones{0};
    long nrhs_max_diff{0};
    long nstep_max_diff{0};
    std::vector<LoomCompStat> comps;
};

// Bit-for-bit equality, except that any NaN matches any NaN: NaN sign and
// payload bits are not compared (user decision 2026-10-18).
inline bool loom_same_bits (double a, double b)
{
    if (a != a && b != b) { return true; }
    std::uint64_t x = 0;
    std::uint64_t y = 0;
    std::memcpy(&x, &a, sizeof x);
    std::memcpy(&y, &b, sizeof y);
    return x == y;
}

inline LoomCompare loom_compare_burns (const amrex::MultiFab& s_ref, const amrex::iMultiFab& n_ref,
                                       const amrex::MultiFab& s_new, const amrex::iMultiFab& n_new,
                                       const plot_t& vars, const amrex::Vector<std::string>& names)
{
    using namespace amrex;
    const int ncomp = s_ref.nComp();
    LoomCompare c;
    c.comps.resize(static_cast<std::size_t>(ncomp));
    for (int n = 0; n < ncomp; ++n) {
        std::string nm = (n < static_cast<int>(names.size())) ? names[n] : ("comp" + std::to_string(n));
        for (char& ch : nm) { if (ch == ' ' || ch == '=' || ch == ',') { ch = '_'; } }
        c.comps[n].name = nm;
    }

    const auto hs = loom_to_host(s_ref);
    const auto hn = loom_to_host(n_ref);
    const auto ls = loom_to_host(s_new);
    const auto ln = loom_to_host(n_new);

    int ib = 0;
    for (MFIter mfi(s_ref); mfi.isValid(); ++mfi, ++ib) {
        const Box bx = mfi.validbox();
        const Dim3 lo = amrex::lbound(bx);
        const Dim3 hi{bx.bigEnd(0) + 1, bx.bigEnd(1) + 1, bx.bigEnd(2) + 1};
        Array4<Real const> const a(hs[ib].data(), lo, hi, ncomp);
        Array4<Real const> const b(ls[ib].data(), lo, hi, ncomp);
        Array4<int const> const an(hn[ib].data(), lo, hi, 2);
        Array4<int const> const bn(ln[ib].data(), lo, hi, 2);
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
            ++c.zones;
            bool same = true;
            for (int n = 0; n < ncomp; ++n) {
                const double x = a(i, j, k, n);
                const double y = b(i, j, k, n);
                if (loom_same_bits(x, y)) { continue; }
                same = false;
                LoomCompStat& st = c.comps[n];
                ++st.differing;
                double ad = std::abs(x - y);
                double rd = (x != 0.0) ? ad / std::abs(x) : std::numeric_limits<double>::infinity();
                if (!std::isfinite(ad)) {   // NaN/inf on either side
                    ad = std::numeric_limits<double>::infinity();
                    rd = std::numeric_limits<double>::infinity();
                }
                if (st.differing == 1 || rd > st.max_rel || (rd == st.max_rel && ad > st.max_abs)) {
                    st.i = i; st.j = j; st.k = k;
                    st.rho = a(i, j, k, vars.irho);
                    st.T = a(i, j, k, vars.itemp);
                }
                st.max_abs = std::max(st.max_abs, ad);
                st.max_rel = std::max(st.max_rel, rd);
            }
            const long dn = std::abs(static_cast<long>(an(i, j, k, 0)) - static_cast<long>(bn(i, j, k, 0)));
            const long ds = std::abs(static_cast<long>(an(i, j, k, 1)) - static_cast<long>(bn(i, j, k, 1)));
            if (dn != 0) { ++c.nrhs_diff_zones; same = false; }
            if (ds != 0) { ++c.nstep_diff_zones; same = false; }
            c.nrhs_max_diff = std::max(c.nrhs_max_diff, dn);
            c.nstep_max_diff = std::max(c.nstep_max_diff, ds);
            if (same) { ++c.identical_zones; }
        });
    }
    return c;
}

// The Loom burn as called from main.cpp (replaces the stock ParReduce when
// unit_test.loom_hsaco is set). Leaves the Loom result in state /
// integrator_n_rhs, stores the Loom failure count in *num_failed_d (the
// device counter main.cpp copies back), and returns the max-n_rhs zone.
inline amrex::ValLocPair<int, burn_t>
loom_burn_main (amrex::MultiFab& state, amrex::iMultiFab& integrator_n_rhs, const plot_t& vars,
                const amrex::Vector<std::string>& names, int* num_failed_d)
{
    using namespace amrex;

    if (ParallelDescriptor::NProcs() > 1) {
        amrex::Abort("Loom launcher: run on one MPI rank");
    }

    const Real dt = unit_test_rp::tmax;
    const int repeats = std::max(1, static_cast<int>(unit_test_rp::loom_repeats));
    LoomBurner lb(unit_test_rp::loom_hsaco, unit_test_rp::loom_kernel, unit_test_rp::loom_ncells);

    int loom_failed = 0;

    if (!unit_test_rp::loom_compare) {
        for (int rep = 0; rep < repeats; ++rep) {
            // the kernel reads only rho, T, ispec_old (unchanged), so
            // repeated burns of the same state are identical
            const LoomRunStats rs = lb.burn(state, integrator_n_rhs, vars, dt);
            amrex::Print() << "LOOM_TIMING run=" << rep << " kind=loom ms=" << rs.ms
                           << " failures=" << rs.num_failed << "\n";
            loom_failed = rs.num_failed;
        }
    } else {
        const int pairs = std::max(2, repeats);
        const BoxArray& ba = state.boxArray();
        const DistributionMapping& dm = state.DistributionMap();
        const int ncomp = state.nComp();

        MultiFab init(ba, dm, ncomp, 0);
        MultiFab::Copy(init, state, 0, 0, ncomp, 0);

        MultiFab s_stock0(ba, dm, ncomp, 0);
        MultiFab s_loom0(ba, dm, ncomp, 0);
        MultiFab s_stock(ba, dm, ncomp, 0);
        MultiFab s_loom(ba, dm, ncomp, 0);
        iMultiFab n_stock0(ba, dm, 2, 0);
        iMultiFab n_loom0(ba, dm, 2, 0);
        iMultiFab n_stock(ba, dm, 2, 0);
        iMultiFab n_loom(ba, dm, 2, 0);

        int stock_failed = 0;
        bool stock_repeat_identical = true;
        bool loom_repeat_identical = true;
        float stock_ms_min = std::numeric_limits<float>::max();
        float loom_ms_min = std::numeric_limits<float>::max();

        for (int p = 0; p < pairs; ++p) {
            MultiFab& ss = (p == 0) ? s_stock0 : s_stock;
            MultiFab& sl = (p == 0) ? s_loom0 : s_loom;
            iMultiFab& ns = (p == 0) ? n_stock0 : n_stock;
            iMultiFab& nl = (p == 0) ? n_loom0 : n_loom;

            MultiFab::Copy(ss, init, 0, 0, ncomp, 0);
            ns.setVal(-1);
            int sf = 0;
            float sms = 0.0f;
            (void) loom_stock_burn(ss, ns, vars, sf, sms);
            amrex::Print() << "LOOM_TIMING run=" << p << " kind=stock ms=" << sms << " failures=" << sf << "\n";
            stock_ms_min = std::min(stock_ms_min, sms);

            MultiFab::Copy(sl, init, 0, 0, ncomp, 0);
            nl.setVal(-1);
            const LoomRunStats rs = lb.burn(sl, nl, vars, dt);
            amrex::Print() << "LOOM_TIMING run=" << p << " kind=loom ms=" << rs.ms
                           << " failures=" << rs.num_failed << "\n";
            loom_ms_min = std::min(loom_ms_min, rs.ms);

            if (p == 0) {
                stock_failed = sf;
                loom_failed = rs.num_failed;
            } else {
                const LoomCompare cs = loom_compare_burns(s_stock0, n_stock0, s_stock, n_stock, vars, names);
                const LoomCompare cl = loom_compare_burns(s_loom0, n_loom0, s_loom, n_loom, vars, names);
                stock_repeat_identical = stock_repeat_identical && cs.identical_zones == cs.zones && sf == stock_failed;
                loom_repeat_identical = loom_repeat_identical && cl.identical_zones == cl.zones &&
                                        rs.num_failed == loom_failed;
            }
        }

        const LoomCompare c = loom_compare_burns(s_stock0, n_stock0, s_loom0, n_loom0, vars, names);
        const bool bit_exact = c.identical_zones == c.zones && stock_failed == loom_failed &&
                               stock_repeat_identical && loom_repeat_identical;

        std::ostringstream line;
        line << std::setprecision(17);
        line << "LOOM_COMPARE zones=" << c.zones << " identical_zones=" << c.identical_zones
             << " nrhs_diff_zones=" << c.nrhs_diff_zones << " nrhs_max_diff=" << c.nrhs_max_diff
             << " nstep_diff_zones=" << c.nstep_diff_zones << " nstep_max_diff=" << c.nstep_max_diff
             << " stock_failures=" << stock_failed << " loom_failures=" << loom_failed
             << " pairs=" << pairs
             << " stock_repeat_identical=" << (stock_repeat_identical ? 1 : 0)
             << " loom_repeat_identical=" << (loom_repeat_identical ? 1 : 0)
             << " stock_ms_min=" << stock_ms_min << " loom_ms_min=" << loom_ms_min
             << " bit_exact=" << (bit_exact ? 1 : 0);
        // per component: name=differing_zones,max_abs,max_rel,i,j,k,rho,T (worst zone by rel)
        for (const LoomCompStat& st : c.comps) {
            line << " " << st.name << "=" << st.differing << "," << st.max_abs << "," << st.max_rel << ","
                 << st.i << "," << st.j << "," << st.k << "," << st.rho << "," << st.T;
        }
        amrex::Print() << line.str() << "\n";

        amrex::Print() << "Loom vs stock, per component (differing zones, max abs, max rel, worst zone i j k rho T):\n";
        for (const LoomCompStat& st : c.comps) {
            if (st.differing == 0) { continue; }
            amrex::Print() << std::setprecision(6) << "  " << st.name << "  " << st.differing << "  "
                           << st.max_abs << "  " << st.max_rel << "  (" << st.i << "," << st.j << "," << st.k
                           << ")  rho=" << st.rho << "  T=" << st.T << "\n";
        }

        if (!bit_exact) {
            amrex::Abort("LOOM_COMPARE: Loom burn is not bit-identical to the stock burn");
        }

        // continue with the Loom result
        MultiFab::Copy(state, s_loom0, 0, 0, ncomp, 0);
        iMultiFab::Copy(integrator_n_rhs, n_loom0, 0, 0, 2, 0);
    }

    Gpu::htod_memcpy(num_failed_d, &loom_failed, sizeof(int));

    return loom_max_nrhs_zone(state, integrator_n_rhs, vars);
}

#endif // MICROPHYSICS_LOOM

#endif
