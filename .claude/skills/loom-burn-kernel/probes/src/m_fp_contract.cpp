#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// Must stay a separate v_mul_f64 + v_add_f64 under --fp-model=precise: the
// correctness references are built with -ffp-contract=off. A v_fma_f64 here
// means Loom contracts, and Phase 6 cannot be bit-close.
namespace probe {
KERNEL void k(double* out, const double* a, const double* b, const double* c) {
  int i = blockIdx.x * blockDim.x + threadIdx.x; out[i] = a[i] * b[i] + c[i]; }
}
