#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// dvstep shape: outer while(true) with retry via continue, inner corrector loop, result code, no return in loops.
namespace probe {
__device__ int step(double* y, int* nrhs) {
  int kflag = 0; int ncf = 0; bool done = false;
  while (true) {
    int m = 0; double acnrm = 0.0; bool conv = false;
    while (true) { y[0] = 0.5 * y[0] + 0.1; ++(*nrhs); ++m; acnrm = y[0] - 0.2; if (acnrm < 1e-3) { conv = true; break; } if (m >= 3) break; }
    if (!conv) { ++ncf; if (ncf >= 10) { kflag = -2; break; } y[0] *= 0.25; continue; }
    done = true; break;
  }
  if (!done && kflag == 0) kflag = -1;
  return kflag;
}
KERNEL void k(double* out, const double* in, int* nrhs) { int i = blockIdx.x * blockDim.x + threadIdx.x; double y[1] = {in[i]}; int n = 0; int kf = step(y, &n); nrhs[i] = n; out[i] = y[0] + kf; }
}
