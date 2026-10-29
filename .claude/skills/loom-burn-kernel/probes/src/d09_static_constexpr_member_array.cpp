#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// static constexpr member array (expected rejected).
namespace probe {
struct T { static constexpr double c[3] = {1.0, 2.0, 3.0}; };
KERNEL void k(double* out, const double* in, const int* sp) { int i = blockIdx.x * blockDim.x + threadIdx.x; out[i] = in[i] * T::c[sp[i]]; }
}
