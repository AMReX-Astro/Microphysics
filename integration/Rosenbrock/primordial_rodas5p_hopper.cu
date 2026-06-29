#include <primordial_rodas5p_hopper.H>

#ifdef AMREX_USE_CUDA

#include <cuda_runtime.h>

#include <string>

#include <AMReX.H>

#include <actual_rhs.H>
#include <burner.H>
#include <integrator_data.H>
#include <linpack.H>
#include <primordial_rodas5p_hopper_generated.H>
#include <rosenbrock_tableau.H>

namespace primordial_rodas5p_hopper {

namespace {

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
T* advance_bytes (unsigned char*& ptr, const std::size_t bytes)
{
    T* out = reinterpret_cast<T*>(ptr);
    ptr += bytes;
    return out;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
bool close_enough (const amrex::Real a, const amrex::Real b)
{
    constexpr amrex::Real rel_tol = 1.0e-13_rt;
    constexpr amrex::Real abs_tol = 1.0e-99_rt;
    return std::abs(a - b) <= abs_tol + rel_tol * amrex::max(std::abs(a), std::abs(b));
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
bool generated_network_matches_scalar (const burn_t& state,
                                       const amrex::Real* rhs_tmp,
                                       const amrex::Real* jac_tmp,
                                       const int local_zone)
{
    RArray1D scalar_rhs{};
    RArray2D scalar_jac{};

    actual_rhs(state, scalar_rhs);
    actual_jac(state, scalar_jac);

    for (int n = 1; n <= neq; ++n) {
        const amrex::Real generated_rhs = rhs_tmp[(n - 1) * zones_per_block + local_zone];
        if (! close_enough(scalar_rhs(n), generated_rhs)) {
            return false;
        }
    }

    for (int j = 1; j <= neq; ++j) {
        for (int i = 1; i <= neq; ++i) {
            const amrex::Real generated_jac =
                jac_tmp[((i - 1) * neq + (j - 1)) * zones_per_block + local_zone];
            if (! close_enough(scalar_jac(i, j), generated_jac)) {
                return false;
            }
        }
    }

    return true;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real& shared_matrix (amrex::Real* matrix, const int i, const int j, const int local_zone)
{
    return matrix[((i - 1) * neq + (j - 1)) * zones_per_block + local_zone];
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real shared_matrix_value (const amrex::Real* matrix, const int i, const int j, const int local_zone)
{
    return matrix[((i - 1) * neq + (j - 1)) * zones_per_block + local_zone];
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void build_rodas5p_matrix (amrex::Real* matrix, const amrex::Real h, const int local_zone)
{
    const amrex::Real fac = 1.0_rt / (h * rosenbrock::rodas5p_tableau::gamma);

    for (int j = 1; j <= neq; ++j) {
        for (int i = 1; i <= neq; ++i) {
            amrex::Real value = -shared_matrix_value(matrix, i, j, local_zone);
            if (i == j) {
                value += fac;
            }
            shared_matrix(matrix, i, j, local_zone) = value;
        }
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
bool rodas5p_matrix_matches_scalar (const burn_t& state,
                                    const amrex::Real* matrix,
                                    const amrex::Real h,
                                    const int local_zone)
{
    RArray2D scalar_jac{};
    actual_jac(state, scalar_jac);

    const amrex::Real fac = 1.0_rt / (h * rosenbrock::rodas5p_tableau::gamma);
    for (int j = 1; j <= neq; ++j) {
        for (int i = 1; i <= neq; ++i) {
            amrex::Real expected = -scalar_jac(i, j);
            if (i == j) {
                expected += fac;
            }

            if (! close_enough(expected, shared_matrix_value(matrix, i, j, local_zone))) {
                return false;
            }
        }
    }

    return true;
}

template <bool allow_pivot>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
bool first_stage_solve_matches_scalar_impl (const burn_t& state,
                                            const amrex::Real* rhs_tmp,
                                            const amrex::Real* matrix,
                                            const amrex::Real h,
                                            const int local_zone)
{
    RArray1D generated_rhs{};
    RArray2D generated_matrix{};
    IArray1D generated_pivot{};

    RArray1D scalar_rhs{};
    RArray2D scalar_matrix{};
    IArray1D scalar_pivot{};

    actual_rhs(state, scalar_rhs);
    actual_jac(state, scalar_matrix);

    const amrex::Real fac = 1.0_rt / (h * rosenbrock::rodas5p_tableau::gamma);

    for (int n = 1; n <= neq; ++n) {
        generated_rhs(n) = rhs_tmp[(n - 1) * zones_per_block + local_zone];
    }

    for (int j = 1; j <= neq; ++j) {
        for (int i = 1; i <= neq; ++i) {
            generated_matrix(i, j) = shared_matrix_value(matrix, i, j, local_zone);

            scalar_matrix(i, j) = -scalar_matrix(i, j);
            if (i == j) {
                scalar_matrix(i, j) += fac;
            }
        }
    }

    int generated_info = 0;
    int scalar_info = 0;
    dgefa<neq, allow_pivot>(generated_matrix, generated_pivot, generated_info);
    dgefa<neq, allow_pivot>(scalar_matrix, scalar_pivot, scalar_info);

    if (generated_info != scalar_info) {
        return false;
    }
    if (generated_info != 0) {
        return true;
    }

    dgesl<neq, allow_pivot>(generated_matrix, generated_pivot, generated_rhs);
    dgesl<neq, allow_pivot>(scalar_matrix, scalar_pivot, scalar_rhs);

    for (int n = 1; n <= neq; ++n) {
        if (! close_enough(scalar_rhs(n), generated_rhs(n))) {
            return false;
        }
    }

    return true;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
bool first_stage_solve_matches_scalar (const burn_t& state,
                                       const amrex::Real* rhs_tmp,
                                       const amrex::Real* matrix,
                                       const amrex::Real h,
                                       const int local_zone)
{
    if (integrator_rp::linalg_do_pivoting == 1) {
        return first_stage_solve_matches_scalar_impl<true>(state, rhs_tmp, matrix, h, local_zone);
    }
    return first_stage_solve_matches_scalar_impl<false>(state, rhs_tmp, matrix, h, local_zone);
}

__global__
void burn_kernel (burn_t* states, const amrex::Real dt, const int nstates)
{
    extern __shared__ unsigned char shared_mem[];

    unsigned char* ptr = shared_mem;
    [[maybe_unused]] amrex::Real* y =
        advance_bytes<amrex::Real>(ptr, bytes_for_reals(neq * zones_per_block));
    [[maybe_unused]] amrex::Real* ynew =
        advance_bytes<amrex::Real>(ptr, bytes_for_reals(neq * zones_per_block));
    [[maybe_unused]] amrex::Real* work =
        advance_bytes<amrex::Real>(ptr, bytes_for_reals(neq * zones_per_block));
    [[maybe_unused]] amrex::Real* ak =
        advance_bytes<amrex::Real>(ptr, bytes_for_reals(8 * neq * zones_per_block));
    [[maybe_unused]] amrex::Real* jac =
        advance_bytes<amrex::Real>(ptr, bytes_for_reals(neq * neq * zones_per_block));
    [[maybe_unused]] amrex::Real* rhs_tmp =
        advance_bytes<amrex::Real>(ptr, bytes_for_reals(neq * zones_per_block));
    [[maybe_unused]] amrex::Real* x_scratch =
        advance_bytes<amrex::Real>(ptr, bytes_for_reals(scratch_count * zones_per_block));
    [[maybe_unused]] int* pivot =
        advance_bytes<int>(ptr, bytes_for_ints(neq * zones_per_block));
    [[maybe_unused]] int* status =
        advance_bytes<int>(ptr, bytes_for_ints(4 * zones_per_block));

    const int warp = threadIdx.x / 32;
    const int lane = threadIdx.x % 32;
    const int local_zone = lane;
    const int zone = blockIdx.x * zones_per_block + local_zone;

    // Prototype milestone 1: establish the Hopper batch launch and comparison
    // harness. The cooperative RHS/Jacobian and shared-memory LU will replace
    // this scalar lane path while preserving the public launcher.
    if (warp == 0 && local_zone < zones_per_block && zone < nstates) {
        primordial_rodas5p_hopper_generated::actual_rhs_scratch(states[zone],
                                                                 rhs_tmp,
                                                                 x_scratch,
                                                                 local_zone,
                                                                 zones_per_block);
        primordial_rodas5p_hopper_generated::actual_jac_scratch(states[zone],
                                                                 jac,
                                                                 x_scratch,
                                                                 local_zone,
                                                                 zones_per_block);
        const bool generated_ok = generated_network_matches_scalar(states[zone],
                                                                    rhs_tmp,
                                                                    jac,
                                                                    local_zone);
        build_rodas5p_matrix(jac, dt, local_zone);
        const bool matrix_ok = rodas5p_matrix_matches_scalar(states[zone],
                                                             jac,
                                                             dt,
                                                             local_zone);
        const bool solve_ok = first_stage_solve_matches_scalar(states[zone],
                                                               rhs_tmp,
                                                               jac,
                                                               dt,
                                                               local_zone);
        status[local_zone] = (generated_ok && matrix_ok && solve_ok) ? 1 : 0;
        burner(states[zone], dt);
        if (! generated_ok || ! matrix_ok || ! solve_ok) {
            states[zone].success = false;
            states[zone].error_code = IERR_BAD_INPUTS;
        }
    }
}

void check_cuda (const cudaError_t err, const char* what)
{
    if (err != cudaSuccess) {
        amrex::Abort(std::string(what) + ": " + cudaGetErrorString(err));
    }
}

} // namespace

void burn_device (burn_t* states, const amrex::Real dt, const int nstates)
{
    if (nstates <= 0) {
        return;
    }

    static bool attributes_set = false;
    if (! attributes_set) {
        check_cuda(cudaFuncSetAttribute(burn_kernel,
                                        cudaFuncAttributeMaxDynamicSharedMemorySize,
                                        static_cast<int>(shared_bytes())),
                   "cudaFuncSetAttribute(MaxDynamicSharedMemorySize)");
        check_cuda(cudaFuncSetAttribute(burn_kernel,
                                        cudaFuncAttributePreferredSharedMemoryCarveout,
                                        cudaSharedmemCarveoutMaxShared),
                   "cudaFuncSetAttribute(PreferredSharedMemoryCarveout)");
        attributes_set = true;
    }

    const int blocks = (nstates + zones_per_block - 1) / zones_per_block;
    burn_kernel<<<blocks, threads_per_block, shared_bytes()>>>(states, dt, nstates);
    check_cuda(cudaGetLastError(), "primordial Rodas5P Hopper kernel launch");
}

} // namespace primordial_rodas5p_hopper

#endif
