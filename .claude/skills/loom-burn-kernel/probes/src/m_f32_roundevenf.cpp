#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
#include <loomcxx/scalar.h>
// OCML BUILTIN_RINT_F32 (cbrt's exponent split) via loom::scalar::roundevenf on f32 (__builtin_rintf does not import).
namespace probe {
KERNEL void k(float* out, const float* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; out[i] = loom::scalar::roundevenf(in[i]); }
}
