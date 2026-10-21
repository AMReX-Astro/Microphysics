#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// static function, double& out-param bound to an automatic local (expected admitted by #1993).
namespace probe {
static void twice(double& a, double& b) { a = 2.0 * a; b = a + 1.0; }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; double a = in[i]; double b = 0.0; twice(a, b); out[i] = a + b; }
}
