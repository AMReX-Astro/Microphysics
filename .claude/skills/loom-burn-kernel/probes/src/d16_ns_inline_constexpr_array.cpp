#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// network_properties-style table as a C++17 inline variable: inline constexpr double aion[N].
namespace probe {
inline constexpr double aion[3] = {4.0, 12.0, 16.0};
KERNEL void k(double* out, const double* in, const int* sp) { int i = blockIdx.x * blockDim.x + threadIdx.x; out[i] = in[i] / aion[sp[i]]; }
}
