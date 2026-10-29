#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
#include <loomcxx/math.h>
// f64 division with the arcp permission (loom::divide_reciprocal, hrx-system #2131).
namespace probe {
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; out[i] = loom::divide_reciprocal(in[2*i], in[2*i+1]); }
}
