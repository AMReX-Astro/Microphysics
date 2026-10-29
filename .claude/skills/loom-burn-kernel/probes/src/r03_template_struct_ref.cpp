#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// template <typename BurnT> void f(BurnT& state) on a local struct (burner/integrator pattern).
namespace probe {
struct Burn { double rho; double T; double xn[13]; int n_rhs; bool success; };
template <typename BurnT> __device__ inline void step(BurnT& s, double dt) { s.T += dt * s.rho; for (int n = 0; n < 13; ++n) s.xn[n] *= 0.5; ++s.n_rhs; s.success = s.T > 0.0; }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; Burn s{}; s.rho = in[i]; s.T = 1.0; step(s, 0.1); out[i] = s.T + s.xn[0] + s.n_rhs; }
}
