#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// reference bound to global memory (expected rejected: needs an automatic SSA object).
namespace probe {
static void twice(double& a) { a = 2.0 * a; }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; out[i] = in[i]; twice(out[i]); }
}
