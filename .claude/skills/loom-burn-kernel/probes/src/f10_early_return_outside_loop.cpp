#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// sneut5: if (temp < 1e7) return 0; before any loop.
namespace probe {
__device__ double sn(double t) { if (t < 1.0e7) { return 0.0; } double s = 0.0; for (int n = 0; n < 4; ++n) s += t * n; return s; }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; out[i] = sn(in[i]); }
}
