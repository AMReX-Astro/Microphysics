#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// double assigned on both branches after declaration (`double r; if (..) r = a; else r = b;`).
namespace probe {
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; double r; if (in[i] > 0.0) r = in[i] * 2.0; else r = -in[i]; out[i] = r; }
}
