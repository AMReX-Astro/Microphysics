#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// double& r = s.x; inside a function (not a parameter).
namespace probe {
struct S { double x; double y; };
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; S s{}; s.x = in[i]; double& r = s.x; r += 1.0; out[i] = s.x; }
}
