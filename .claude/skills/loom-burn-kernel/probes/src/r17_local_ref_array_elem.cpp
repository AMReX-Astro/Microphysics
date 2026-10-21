#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// local reference bound to a local array element (`double& r = a[2]`).
namespace probe {
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; double a[4] = {}; double& r = a[2]; r = in[i]; out[i] = a[2]; }
}
