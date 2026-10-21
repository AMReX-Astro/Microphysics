#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// local array without initializer (`double a[14];`), filled by a loop.
namespace probe {
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; double a[14]; for (int n = 0; n < 14; ++n) a[n] = in[i] * n; double s = 0.0; for (int n = 0; n < 14; ++n) s += a[n]; out[i] = s; }
}
