#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// block-scope static constexpr array inside a device function, runtime index.
namespace probe {
__device__ double lookup(int n) { static constexpr double t[3] = {4.0, 12.0, 16.0}; return t[n]; }
KERNEL void k(double* out, const double* in, const int* sp) { int i = blockIdx.x * blockDim.x + threadIdx.x; out[i] = in[i] / lookup(sp[i] % 3); }
}
