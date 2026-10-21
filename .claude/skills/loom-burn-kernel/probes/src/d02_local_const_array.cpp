#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// ifermi12 a1/a2/b1 and the aion replacement: local const array, runtime index.
namespace probe {
KERNEL void k(double* out, const double* in, const int* sp) { int i = blockIdx.x * blockDim.x + threadIdx.x; const double aion[13] = {4.0, 12.0, 16.0, 20.0, 24.0, 28.0, 32.0, 36.0, 40.0, 44.0, 48.0, 52.0, 56.0}; out[i] = in[i] / aion[sp[i]]; }
}
