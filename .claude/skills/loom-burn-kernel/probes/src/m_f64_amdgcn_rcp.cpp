#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// OCML BUILTIN_AMDGPU_RCP_F64 (MATH_FAST_RCP / MATH_FAST_DIV seed): __builtin_amdgcn_rcp on f64.
namespace probe {
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; out[i] = __builtin_amdgcn_rcp(in[i]); }
}
