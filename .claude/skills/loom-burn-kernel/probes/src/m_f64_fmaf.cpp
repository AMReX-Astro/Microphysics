#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
#include <loomcxx/scalar.h>
namespace probe {
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; out[i] = loom::scalar::fmaf(in[3*i], in[3*i+1], in[3*i+2]); }
}
