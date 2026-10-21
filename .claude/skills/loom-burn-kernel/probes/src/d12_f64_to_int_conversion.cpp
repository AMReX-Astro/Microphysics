#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// helmholtz table hash: int((log10(T) - tlo) * tstpi) + 1 (f64 -> i32 conversion).
namespace probe {
constexpr double tlo = 3.0; constexpr double tstpi = 1.0 / 0.05;
KERNEL void k(int* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; out[i] = int((in[i] - tlo) * tstpi) + 1; }
}
