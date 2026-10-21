#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
#include <loomcxx/scalar.h>
// f64 frexp through loom::scalar::frexpf (hrx-system #2131): fraction and i32 exponent.
namespace probe {
KERNEL void k(double* out, const double* in, int* ex) { int i = blockIdx.x * blockDim.x + threadIdx.x; auto r = loom::scalar::frexpf(in[i]); out[i] = r.fraction; ex[i] = r.exponent; }
}
