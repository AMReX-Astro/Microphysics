#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// integrator_wrapper retry: burn_t old_state{state}; ... state = old_state;
namespace probe {
struct Burn { double rho; double T; double xn[13]; double ydot_a[15]; int n_rhs; bool success; };
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; Burn s{}; s.T = in[i]; Burn old{s}; s.T = -1.0; s.success = false; if (!s.success) { s = old; } out[i] = s.T; }
}
