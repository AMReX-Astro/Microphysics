#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// scalar fields mutated by a member function called through a struct reference (`s.bump(dt)`).
namespace probe {
struct Burn { double rho; double T; int n_rhs; __device__ void bump(double dt) { T += dt * rho; ++n_rhs; } };
__device__ static void step(Burn& s, double dt) { s.bump(dt); }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; Burn s{}; s.rho = in[i]; s.T = 1.0; step(s, 0.1); out[i] = s.T + s.n_rhs; }
}
