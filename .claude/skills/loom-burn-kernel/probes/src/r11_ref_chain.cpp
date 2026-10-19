#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// reference parameter forwarded through two static calls (burner -> integrator -> setup).
namespace probe {
struct Burn { double T; double e; int n; };
static void inner(Burn& s) { s.e = 2.0 * s.T; ++s.n; }
static void outer(Burn& s, double dt) { s.T += dt; inner(s); }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; Burn s{}; s.T = in[i]; outer(s, 0.5); out[i] = s.e + s.n; }
}
