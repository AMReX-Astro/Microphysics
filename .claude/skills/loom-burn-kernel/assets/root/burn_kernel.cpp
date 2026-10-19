// Loom kernel root for a helmholtz + VODE burn (see README.md):
// test_react's per-zone burn (unit_test/test_react/react_zones.H do_react),
// one lane per zone of one box.
//
// Root: mpl::burn_zones (C++) -> mpl.burn_zones (Loom), workgroup size 128.
// The workgroup count is fixed at loom-compile time:
//   workgroup_count.x = ceil(nx * ny * nz / 128).
//
// Arguments (pointers and int/long/double scalars only):
//   params        BurnParams (scalars) filled on the host (loom_params.H)
//   f ... dd2i_sav the 14 helmholtz tables (EosTables fields, in order);
//                 EosTables is rebuilt BY VALUE here (a pointer to a struct
//                 of pointers does not import)
//   state         the state box's data pointer (Array4::p of the MultiFab
//                 box), with its jstride/kstride/nstride and the box lower
//                 corner lo_i/lo_j/lo_k and extents nx/ny/nz
//   irho ... irho_hnuc  plot_t component indices (variables.H)
//   n_rhs         the integrator_n_rhs iMultiFab box data (2 components:
//                 n_rhs, n_step), with its own strides; same box as state
//   dt            unit_test.tmax
//   zone_success  one int per zone (linear zone index), written 1/0; the
//                 host reduces it into num_failed
//
// Zone index: idx = blockIdx.x * blockDim.x + threadIdx.x, decomposed as an
// AMReX loop over the box (i fastest, then j, then k). Lanes with
// idx >= nx * ny * nz do nothing.
//
// On the host (no __loom__), burn_zones is an ordinary inline function with
// one extra trailing parameter, the zone index, so the native validator
// (validate_kernel.H) runs exactly this body zone by zone. With hipcc and
// -DMPL_HIP_STANDALONE (assets/standalone/main.hip) it is a HIP __global__
// kernel with the same signature as the Loom root, launched with
// blockDim.x = 128 and gridDim.x = ceil(nx * ny * nz / 128).
#include "burner.H"

#if defined(MPL_GPU_KERNEL)
#define MPL_BURN_ZONE_INDEX_PARAM
#else
#define MPL_BURN_ZONE_INDEX_PARAM , int host_zone_index
#endif

namespace mpl {

// from unit_test/test_react/react_zones.H:13 do_react (NAUX_NET = 0,
// no NONAKA_PLOT) and the main.cpp ParReduce body (failure count -> per-zone
// success flag)
MPL_KERNEL void burn_zones(const BurnParams* params,
                           const Real* f, const Real* dpdf, const Real* ef, const Real* xf,
                           const Real* d, const Real* t, const Real* dt_sav, const Real* dt2_sav,
                           const Real* dti_sav, const Real* dt2i_sav, const Real* dd_sav,
                           const Real* dd2_sav, const Real* ddi_sav, const Real* dd2i_sav,
                           Real* state, long jstride, long kstride, long nstride,
                           int lo_i, int lo_j, int lo_k, int nx, int ny, int nz,
                           int irho, int itemp, int ispec_old, int ispec, int irodot, int irho_hnuc,
                           int* n_rhs, long n_rhs_jstride, long n_rhs_kstride, long n_rhs_nstride,
                           Real dt, int* zone_success MPL_BURN_ZONE_INDEX_PARAM)
{
    EosTables tab{f, dpdf, ef, xf, d, t, dt_sav, dt2_sav, dti_sav, dt2i_sav,
                  dd_sav, dd2_sav, ddi_sav, dd2i_sav};

#if defined(MPL_GPU_KERNEL)
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
#else
    const int idx = host_zone_index;
#endif
    const int ncells = nx * ny * nz;

    if (idx < ncells) {
        const int nxy = nx * ny;
        const int ii = idx % nx;
        const int jj = (idx / nx) % ny;
        const int kk = idx / nxy;

        // Array4 offsets of (i, j, k) relative to the box lower corner
        const long cell = static_cast<long>(ii) + static_cast<long>(jj) * jstride +
                          static_cast<long>(kk) * kstride;
        const long cell_n = static_cast<long>(ii) + static_cast<long>(jj) * n_rhs_jstride +
                            static_cast<long>(kk) * n_rhs_kstride;

        burn_t burn_state{};

        burn_state.rho = state[cell + static_cast<long>(irho) * nstride];
        burn_state.T = state[cell + static_cast<long>(itemp) * nstride];
        for (int n = 0; n < NumSpec; ++n) {
            burn_state.xn[n] = state[cell + static_cast<long>(ispec_old + n) * nstride];
        }

        const IntegratorParams ip = params->integ;
        normalize_abundances_burn(&burn_state, &ip);

        // the integrator doesn't actually care about the initial internal
        // energy.
        burn_state.e = 0.0;

        burn_state.i = lo_i + ii;
        burn_state.j = lo_j + jj;
        burn_state.k = lo_k + kk;

        burn_state.T_fixed = -1.0;

        burner(&burn_state, dt, params, tab);

        for (int n = 0; n < NumSpec; ++n) {
            state[cell + static_cast<long>(ispec + n) * nstride] = burn_state.xn[n];
            state[cell + static_cast<long>(irodot + n) * nstride] =
                (burn_state.xn[n] - state[cell + static_cast<long>(ispec_old + n) * nstride]) / dt;
        }

        state[cell + static_cast<long>(irho_hnuc) * nstride] =
            state[cell + static_cast<long>(irho) * nstride] * burn_state.e / dt;

        n_rhs[cell_n] = burn_state.n_rhs;
        n_rhs[cell_n + n_rhs_nstride] = burn_state.n_step;

        zone_success[idx] = burn_state.success ? 1 : 0;
    }
}

} // namespace mpl

#undef MPL_BURN_ZONE_INDEX_PARAM
