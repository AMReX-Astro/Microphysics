#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// scalar field written through a whole-struct reference (`s.T = ...` in f(Burn&)).
namespace probe {
struct Burn { double rho; double T; double xn[13]; int n_rhs; bool success; };
__device__ static void step(Burn& s, double dt) { s.T = s.T + dt * s.rho; }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; Burn s{}; s.rho = in[i]; s.T = 1.0; step(s, 0.1); out[i] = s.T + s.xn[0] + s.n_rhs; }
}
