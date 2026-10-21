#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// state and a member of state passed as two mutable refs (expected rejected).
namespace probe {
struct S { double e; double T; };
static void upd(S& s, double& t) { s.e += t; t = s.e; }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; S s{}; s.e = in[i]; upd(s, s.T); out[i] = s.e + s.T; }
}
