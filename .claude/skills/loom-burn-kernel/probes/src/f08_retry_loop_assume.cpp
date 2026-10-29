#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// same as f07 with __builtin_assume bounds on the counter.
namespace probe {
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; double h = in[i]; int ntry = 0;
  while (true) { __builtin_assume(0 <= ntry && ntry < 10); ++ntry; double err = h * h; if (err > 1.0 && ntry < 10) { h *= 0.5; continue; } break; } out[i] = h + ntry; }
}
