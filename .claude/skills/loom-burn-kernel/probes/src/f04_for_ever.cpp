#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// for (;;) (expected rejected).
namespace probe {
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; double y = in[i]; int m = 0; for (;;) { y *= 0.5; ++m; if (y < 1.0 || m > 50) break; } out[i] = y; }
}
