#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// const double (&y)[13] (expected rejected: array reference parameters require addressable storage).
namespace probe {
static double first(const double (&y)[13]) { return y[0] + y[12]; }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; double y[13] = {}; y[0] = in[i]; out[i] = first(y); }
}
