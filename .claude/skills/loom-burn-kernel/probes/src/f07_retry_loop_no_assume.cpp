#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// stage/retry loop with break/continue and a counter, no __builtin_assume (SUBRANGE/026 risk).
namespace probe {
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; double h = in[i]; int ntry = 0;
  while (true) { ++ntry; double err = h * h; if (err > 1.0 && ntry < 10) { h *= 0.5; continue; } break; } out[i] = h + ntry; }
}
