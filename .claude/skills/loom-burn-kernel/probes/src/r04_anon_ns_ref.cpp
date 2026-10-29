#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// non-static function in an anonymous namespace (internal linkage) with a reference param.
namespace probe {
namespace { __device__ void twice(double& a) { a = 2.0 * a; } }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; double a = in[i]; twice(a); out[i] = a; }
}
