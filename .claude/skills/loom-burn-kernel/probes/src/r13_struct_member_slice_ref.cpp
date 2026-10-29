#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// pass a member of a local struct as a mutable reference (non-overlapping member slice).
namespace probe {
struct S { double e; double T; };
static void bump(double& t) { t += 1.0; }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; S s{}; s.T = in[i]; bump(s.T); out[i] = s.T + s.e; }
}
