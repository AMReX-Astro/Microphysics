#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// conditionless `for (int m = 0;; ++m)` with continue and break.
namespace probe {
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; double y = in[i]; int m = 0; for (m = 0;; ++m) { if (m < 2) continue; y *= 0.5; if (y < 1.0 || m > 50) break; } out[i] = y + m; }
}
