#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// actual_log_screen shape: early return inside an if nested in a runtime if, then a fallthrough.
namespace probe {
__device__ double f(double x, int skip, double thr) { double r = 0.0; if (skip) { r = x * 0.5; if (r <= thr) { return r; } } r = x * x; return r; }
KERNEL void k(double* out, const double* in, const int* sp) { int i = blockIdx.x * blockDim.x + threadIdx.x; out[i] = f(in[i], sp[i], 1.0); }
}
