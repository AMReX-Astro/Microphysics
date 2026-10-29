#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// array-member elements written through a whole-struct reference (`s.xn[n] = ...`).
namespace probe {
struct Burn { double rho; double T; double xn[13]; int n_rhs; bool success; };
__device__ static void step(Burn& s, double dt) { for (int n = 0; n < 13; ++n) s.xn[n] = dt * n; }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; Burn s{}; s.rho = in[i]; s.T = 1.0; step(s, 0.1); out[i] = s.T + s.xn[0] + s.n_rhs; }
}
