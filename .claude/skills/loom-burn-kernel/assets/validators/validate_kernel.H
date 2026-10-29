// Native bit-exact validation of the kernel root burn_kernel.cpp
// (mpl::burn_zones, host path) against test_react's stock per-zone burn
// (react_zones.H do_react). Host only; included by loom_validate.cpp under
// -DLOOM_VALIDATE_KERNEL. See README.md.
//
// A state buffer is laid out exactly like one Array4 box of test_react's
// state MultiFab (plot_t components) and filled exactly as main.cpp does
// (rho/T log grids, setup_composition/get_xn, ispec_old = max(X, 1e-10)).
// One copy is burned by the stock do_react (through an amrex::Array4 over the
// buffer), the other by mpl::burn_zones, zone by zone; every component of
// every zone, n_rhs/n_step and the success flag are compared bit for bit.
// Boxes: n_cell = 16 and 32 (one box, lo = 0, as main.cpp with
// max_grid_size = n_cell), a sub-box with a nonzero, non-cubic lower corner
// (index decomposition), and n_cell = 16 with the prad limiter on.
// integrator.use_jacobian_caching is forced to 0 (the HIP path) while it
// runs. Sensitivity: dt + 1 ulp on the Loom side must produce mismatches.
#ifndef VALIDATE_KERNEL_H
#define VALIDATE_KERNEL_H

#include <cmath>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

#include <AMReX_Array4.H>
#include <extern_parameters.H>
#include <variables.H>
#include <react_util.H>
#include <react_zones.H>     // stock do_react (the kernel dir is last on -I)
#include <loom_validate.H>
#include <loom_params.H>

#include "burn_kernel.cpp"   // kernel root (this directory)

namespace loomval_kernel {

using amrex::Real;

struct KBox {
    int lo[3];
    int n[3];
};

// main.cpp initialization of one box (global grid of n_cell^3 zones)
inline void init_box(const KBox& b, int n_cell, const plot_t& vars, amrex::Array4<Real> const& s)
{
    Real dlogrho = 0.0_rt;
    Real dlogT = 0.0_rt;
    if (n_cell > 1) {
        dlogrho = (std::log10(unit_test_rp::dens_max) - std::log10(unit_test_rp::dens_min)) / static_cast<amrex::Real>(n_cell - 1);
        dlogT   = (std::log10(unit_test_rp::temp_max) - std::log10(unit_test_rp::temp_min)) / static_cast<amrex::Real>(n_cell - 1);
    }
    init_t comp_data = setup_composition(n_cell);
    for (int k = b.lo[2]; k < b.lo[2] + b.n[2]; ++k) {
        for (int j = b.lo[1]; j < b.lo[1] + b.n[1]; ++j) {
            for (int i = b.lo[0]; i < b.lo[0] + b.n[0]; ++i) {
                s(i, j, k, vars.itemp) =
                    std::pow(10.0_rt, (std::log10(unit_test_rp::temp_min) + static_cast<amrex::Real>(j)*dlogT));
                s(i, j, k, vars.irho) =
                    std::pow(10.0_rt, (std::log10(unit_test_rp::dens_min) + static_cast<amrex::Real>(i)*dlogrho));
                Real xn[NumSpec];
                get_xn(k, comp_data, xn, unit_test_rp::uniform_xn);
                for (int n = 0; n < NumSpec; n++) {
                    s(i, j, k, vars.ispec_old+n) = amrex::max(xn[n], 1.e-10_rt);
                }
            }
        }
    }
}

struct RunStats {
    long zones{0};
    long stock_fail{0};
    long loom_fail{0};
};

// Burn one box both ways and compare. dt_perturb != 0 shifts the Loom dt by
// one ulp (sensitivity check).
inline RunStats run_box(LoomCheck& c, const std::string& tag, const KBox& b, int n_cell,
                        const plot_t& vars, int dt_perturb)
{
    const amrex::Dim3 lo{b.lo[0], b.lo[1], b.lo[2]};
    const amrex::Dim3 hi{b.lo[0] + b.n[0], b.lo[1] + b.n[1], b.lo[2] + b.n[2]};  // exclusive end
    const long ncell_box = static_cast<long>(b.n[0]) * b.n[1] * b.n[2];
    const int ncomp = vars.n_plot_comps;

    std::vector<Real> buf_stock(static_cast<std::size_t>(ncell_box * ncomp), 0.0_rt);
    std::vector<int> nrhs_stock(static_cast<std::size_t>(ncell_box * 2), -7);
    amrex::Array4<Real> s_stock(buf_stock.data(), lo, hi, ncomp);
    amrex::Array4<int> n_stock(nrhs_stock.data(), lo, hi, 2);
    init_box(b, n_cell, vars, s_stock);

    std::vector<Real> buf_loom = buf_stock;
    std::vector<int> nrhs_loom = nrhs_stock;
    amrex::Array4<Real> s_loom(buf_loom.data(), lo, hi, ncomp);
    amrex::Array4<int> n_loom(nrhs_loom.data(), lo, hi, 2);
    std::vector<int> success_loom(static_cast<std::size_t>(ncell_box), -1);
    std::vector<int> success_stock(static_cast<std::size_t>(ncell_box), -1);

    // stock: AMReX loop order (i fastest)
    long lin = 0;
    for (int k = lo.z; k < hi.z; ++k) {
        for (int j = lo.y; j < hi.y; ++j) {
            for (int i = lo.x; i < hi.x; ++i) {
                burn_t burn_state;
                bool success = do_react(i, j, k, s_stock, burn_state, n_stock, vars);
                success_stock[lin++] = success ? 1 : 0;
            }
        }
    }

    // Loom root, host path
    mpl::BurnParams bp{};
    fill_burn_params(&bp);
    mpl::EosTables tab{};
    fill_eos_tables(&tab);
    Real dt = unit_test_rp::tmax;
    if (dt_perturb != 0) { dt = std::nextafter(dt, 1.0e300); }
    const long jstride = s_loom.get_stride<1>();
    const long kstride = s_loom.get_stride<2>();
    const long nstride = s_loom.get_stride<3>();
    const long n_jstride = n_loom.get_stride<1>();
    const long n_kstride = n_loom.get_stride<2>();
    const long n_nstride = n_loom.get_stride<3>();
    for (int idx = 0; idx < static_cast<int>(ncell_box); ++idx) {
        mpl::burn_zones(&bp,
                        tab.f, tab.dpdf, tab.ef, tab.xf, tab.d, tab.t, tab.dt_sav, tab.dt2_sav,
                        tab.dti_sav, tab.dt2i_sav, tab.dd_sav, tab.dd2_sav, tab.ddi_sav, tab.dd2i_sav,
                        s_loom.p, jstride, kstride, nstride,
                        b.lo[0], b.lo[1], b.lo[2], b.n[0], b.n[1], b.n[2],
                        vars.irho, vars.itemp, vars.ispec_old, vars.ispec, vars.irodot, vars.irho_hnuc,
                        n_loom.p, n_jstride, n_kstride, n_nstride,
                        dt, success_loom.data(), idx);
    }
    // a lane past the end must not touch anything
    mpl::burn_zones(&bp,
                    tab.f, tab.dpdf, tab.ef, tab.xf, tab.d, tab.t, tab.dt_sav, tab.dt2_sav,
                    tab.dti_sav, tab.dt2i_sav, tab.dd_sav, tab.dd2_sav, tab.ddi_sav, tab.dd2i_sav,
                    s_loom.p, jstride, kstride, nstride,
                    b.lo[0], b.lo[1], b.lo[2], b.n[0], b.n[1], b.n[2],
                    vars.irho, vars.itemp, vars.ispec_old, vars.ispec, vars.irodot, vars.irho_hnuc,
                    n_loom.p, n_jstride, n_kstride, n_nstride,
                    dt, success_loom.data(), static_cast<int>(ncell_box));

    RunStats st;
    lin = 0;
    for (int k = lo.z; k < hi.z; ++k) {
        for (int j = lo.y; j < hi.y; ++j) {
            for (int i = lo.x; i < hi.x; ++i) {
                const std::string z = tag + "(" + std::to_string(i) + "," + std::to_string(j) + "," +
                                      std::to_string(k) + ")";
                for (int n = 0; n < ncomp; ++n) {
                    loom_check_bits(c, (z + ".comp" + std::to_string(n)).c_str(),
                                    s_stock(i, j, k, n), s_loom(i, j, k, n));
                }
                loom_check_int(c, (z + ".n_rhs").c_str(), n_stock(i, j, k, 0), n_loom(i, j, k, 0));
                loom_check_int(c, (z + ".n_step").c_str(), n_stock(i, j, k, 1), n_loom(i, j, k, 1));
                loom_check_int(c, (z + ".success").c_str(), success_stock[lin], success_loom[lin]);
                ++st.zones;
                if (success_stock[lin] == 0) { ++st.stock_fail; }
                if (success_loom[lin] == 0) { ++st.loom_fail; }
                ++lin;
            }
        }
    }
    return st;
}

inline void print_stats(const std::string& what, const RunStats& st)
{
    amrex::Print() << "LOOM_VALIDATE kernel " << what << " zones=" << st.zones
                   << " stock_failures=" << st.stock_fail << " loom_failures=" << st.loom_fail << "\n";
}

inline long validate_kernel()
{
    LoomCheck c{"kernel"};

    const bool saved_caching = integrator_rp::use_jacobian_caching;
    const Real saved_rho_c = eos_rp::prad_limiter_rho_c;
    const Real saved_delta_rho = eos_rp::prad_limiter_delta_rho;
    integrator_rp::use_jacobian_caching = false;   // the HIP (no caching) path

    const plot_t vars = init_variables();

    for (int n_cell : {16, 32}) {
        const KBox b{{0, 0, 0}, {n_cell, n_cell, n_cell}};
        const RunStats st = run_box(c, "box" + std::to_string(n_cell), b, n_cell, vars, 0);
        print_stats("n_cell=" + std::to_string(n_cell) + " box=0..n_cell-1", st);
    }

    {
        // sub-box of the n_cell = 16 grid: nonzero, non-cubic lower corner
        const KBox b{{3, 6, 9}, {11, 7, 5}};
        const RunStats st = run_box(c, "subbox16", b, 16, vars, 0);
        print_stats("n_cell=16 subbox lo=(3,6,9) n=(11,7,5)", st);
    }

    {
        // prad limiter on (both parameters must be > 0 to enable it)
        eos_rp::prad_limiter_rho_c = 1.0e7_rt;
        eos_rp::prad_limiter_delta_rho = 3.0e6_rt;
        const KBox b{{0, 0, 0}, {16, 16, 16}};
        const RunStats st = run_box(c, "prad16", b, 16, vars, 0);
        print_stats("n_cell=16 prad_limiter_rho_c=1e7 prad_limiter_delta_rho=3e6", st);
        eos_rp::prad_limiter_rho_c = saved_rho_c;
        eos_rp::prad_limiter_delta_rho = saved_delta_rho;
    }

    {
        // sensitivity: Loom dt + 1 ulp must be detected. Use the full grid:
        // a corner (the coldest, least dense zones) may not burn at all for
        // a gentle network (triple_alpha_plus_cago on inputs_3alpha), and then
        // a 1-ulp dt change leaves every output bit unchanged.
        LoomCheck cs{"kernel_sensitivity"};
        cs.reported = 20;  // expected mismatches: do not print them
        const KBox b{{0, 0, 0}, {16, 16, 16}};
        run_box(cs, "sens", b, 16, vars, 1);
        amrex::Print() << "LOOM_VALIDATE kernel_sensitivity (dt +1ulp, expected nonzero) checks="
                       << cs.checks << " mismatches=" << cs.mismatches << "\n";
        if (cs.mismatches == 0) {
            amrex::Print() << "LOOM_MISMATCH kernel_sensitivity: perturbation went undetected\n";
            c.mismatches += 1;
        }
    }

    integrator_rp::use_jacobian_caching = saved_caching;
    eos_rp::prad_limiter_rho_c = saved_rho_c;
    eos_rp::prad_limiter_delta_rho = saved_delta_rho;
    return loom_report(c);
}

} // namespace loomval_kernel

using loomval_kernel::validate_kernel;

#endif
