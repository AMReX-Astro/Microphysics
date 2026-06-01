#include <primordial_rodas5p_hopper.H>

#ifdef AMREX_USE_CUDA

#include <cuda_runtime.h>

#include <string>

#include <AMReX.H>

#include <actual_rhs.H>
#include <burner.H>
#include <integrator_data.H>
#include <primordial_rodas5p_hopper_generated.H>

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
        status[local_zone] = generated_ok ? 1 : 0;
        burner(states[zone], dt);
        if (! generated_ok) {
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
