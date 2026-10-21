#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
#include <loomcxx/scalar.h>
// f64 division with the afn permission (loom::scalar::approximate::divf).
namespace probe {
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; out[i] = loom::scalar::approximate::divf(in[2*i], in[2*i+1]); }
}
