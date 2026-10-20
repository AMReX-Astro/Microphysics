#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// retry loop nested in a step loop, each with several breaks (dvstep corrector shape)
namespace probe {
struct St { double t; double h; double y[4]; int n_step; int nfail; int result; };
__device__ static int step(St* s) { double e = s->h * s->y[0]; if (e > 2.0) { s->nfail += 1; s->h *= 0.25; return -1; } for (int n = 0; n < 4; ++n) s->y[n] *= 0.9; s->t += s->h; return 0; }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; St s{}; s.h = in[i]; for (int n = 0; n < 4; ++n) s.y[n] = in[i] + n;
  while (true) {
    int m = 0; double del = 1.0;
    while (true) { del *= s.y[m % 4] * 0.3; m += 1; if (del < 1.0e-3) break; if (m >= 3) { s.nfail += 1; break; } if (del > 1.0e3) { s.nfail += 2; break; } }
    if (s.nfail >= 10) { s.result = -1; break; }
    if (m >= 3) { s.h *= 0.5; continue; }
    s.t += s.h; s.n_step += 1;
    if (s.t >= 10.0 || s.n_step >= 200) { s.result = 1; break; }
  }
  out[i] = s.t + s.result; }
}
