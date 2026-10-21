#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// f24 with the two conditions merged into one if: one early return, then a fallthrough.
namespace probe {
__device__ double f(double x, int skip, double thr) { double r = x * 0.5; if (skip && r <= thr) { return r; } r = x * x; return r; }
KERNEL void k(double* out, const double* in, const int* sp) { int i = blockIdx.x * blockDim.x + threadIdx.x; out[i] = f(in[i], sp[i], 1.0); }
}
