#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// OCML spelling: __builtin_frexp with an int out-parameter, then __builtin_ldexp.
namespace probe {
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; int e; double m = __builtin_frexp(in[i], &e); out[i] = __builtin_ldexp(m, e - 3); }
}
