#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
#include <loomcxx/scalar.h>
// f64 ldexp through loom::scalar::ldexpf (hrx-system #2131).
namespace probe {
KERNEL void k(double* out, const double* in, const int* n) { int i = blockIdx.x * blockDim.x + threadIdx.x; out[i] = loom::scalar::ldexpf(in[i], n[i]); }
}
