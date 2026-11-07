#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// r08 with a real alias: the callee writes s.T through t, then reads s.T through s.
// Copy-in/copy-out of t gives s.e = in instead of 4*in.
namespace probe {
struct S { double e; double T; };
static void upd(S& s, double& t) { t = 3.0 * s.e; s.e += s.T; }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; S s{}; s.e = in[i]; upd(s, s.T); out[i] = s.e + 10.0 * s.T; }
}
