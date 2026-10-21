#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// Nordsieck advance: for k = NQ..1, for j = k..NQ with runtime NQ (vstate.NQ).
namespace probe {
KERNEL void k(double* out, const double* in, const int* nq_in) { int i = blockIdx.x * blockDim.x + threadIdx.x;
  double yh[14 * 6] = {}; for (int n = 0; n < 14 * 6; ++n) yh[n] = in[i] + n; int NQ = nq_in[i];
  for (int kk = NQ; kk >= 1; --kk) { for (int j = kk; j <= NQ; ++j) { for (int e = 0; e < 14; ++e) { yh[e * 6 + (j - 1)] += yh[e * 6 + j]; } } }
  out[i] = yh[0] + yh[13 * 6 + 5]; }
}
