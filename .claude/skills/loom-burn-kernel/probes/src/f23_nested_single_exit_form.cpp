#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// x2 rewritten in the single-exit form (done flags, guarded tails), the current rewrite
namespace probe {
struct St { double t; double h; double y[4]; int n_step; int nfail; int result; };
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; St s{}; s.h = in[i]; for (int n = 0; n < 4; ++n) s.y[n] = in[i] + n;
  bool done = false;
  while (!done) {
    int m = 0; double del = 1.0; bool idone = false;
    while (!idone) { del *= s.y[m % 4] * 0.3; m += 1; if (del < 1.0e-3) { idone = true; } else if (m >= 3) { s.nfail += 1; idone = true; } else if (del > 1.0e3) { s.nfail += 2; idone = true; } }
    if (s.nfail >= 10) { s.result = -1; done = true; }
    else if (m >= 3) { s.h *= 0.5; }
    else { s.t += s.h; s.n_step += 1; if (s.t >= 10.0 || s.n_step >= 200) { s.result = 1; done = true; } }
  }
  out[i] = s.t + s.result; }
}
