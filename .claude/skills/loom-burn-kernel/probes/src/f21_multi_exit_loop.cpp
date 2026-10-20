#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// dvode step-loop shape: while(true), 4 breaks + 1 continue inside data-dependent branches, side effects before each exit
namespace probe {
struct St { double t; double h; double y[4]; int n_step; int nfail; int result; };
__device__ static int step(St* s) { double e = s->h * s->y[0]; if (e > 2.0) { s->nfail += 1; s->h *= 0.25; return -1; } for (int n = 0; n < 4; ++n) s->y[n] *= 0.9; s->t += s->h; return 0; }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; St s{}; s.h = in[i]; for (int n = 0; n < 4; ++n) s.y[n] = in[i] + n; bool first = true;
  while (true) {
    if (!first) { if (s.n_step >= 100) { s.result = -1; for (int n = 0; n < 4; ++n) s.y[n] = 0.0; break; } if (s.h < 1.0e-12) { s.result = -2; break; } }
    first = false;
    int kflag = step(&s);
    if (kflag < 0) { if (s.nfail >= 7) { s.result = -3; break; } continue; }
    s.n_step += 1;
    if (s.t >= 10.0) { s.result = 1; break; }
    s.h *= 1.5;
  }
  out[i] = s.t + s.result + s.y[0]; }
}
