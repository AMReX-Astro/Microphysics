#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// const Burn& read-only, reading a member array (eos/rhs read pattern).
namespace probe {
struct Burn { double rho; double xn[13]; };
static double sum(const Burn& s) { double t = 0.0; for (int n = 0; n < 13; ++n) t += s.xn[n]; return t * s.rho; }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; Burn s{}; s.rho = in[i]; s.xn[3] = 1.0; out[i] = sum(s); }
}
