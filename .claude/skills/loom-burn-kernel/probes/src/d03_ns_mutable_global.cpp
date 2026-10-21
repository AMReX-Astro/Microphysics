#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// AMREX_GPU_MANAGED runtime parameter as a namespace-scope global (expected rejected).
namespace probe {
namespace integrator_rp { double rtol_spec = 1.0e-12; }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; out[i] = in[i] * integrator_rp::rtol_spec; }
}
