#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// AMREX_INLINE / AMREX_GPU_HOST_DEVICE style: inline with external linkage, not static.
namespace probe {
__device__ inline void twice(double& a) { a = 2.0 * a; }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; double a = in[i]; twice(a); out[i] = a; }
}
