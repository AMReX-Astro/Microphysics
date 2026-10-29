#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// helmholtz actual_eos: for (iter = 1..100) Newton loop with break on convergence.
namespace probe {
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; double T = 1.0; const double want = in[i]; bool converged = false;
  for (int iter = 1; iter <= 100; ++iter) { double f = T * T * T - want; double df = 3.0 * T * T; if (converged) break; double dT = f / df; T -= dT; if (dT < 1e-12 * T && -dT < 1e-12 * T) converged = true; } out[i] = T; }
}
