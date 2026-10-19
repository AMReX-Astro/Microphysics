#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// user-declared constructor with arguments (autodiff Dual(x), rate_t{} etc.).
namespace probe {
struct D { double v; double d; __device__ D(double x) : v(x), d(0.0) {} };
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; D x(in[i]); out[i] = x.v + x.d; }
}
