#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// namespace-scope constexpr arrays read with a constant-bound loop counter (sum of X/A, Ye).
namespace probe {
constexpr double aion[4] = {1.0, 4.0, 12.0, 16.0};
constexpr double zion[4] = {1.0, 2.0, 6.0, 8.0};
KERNEL void k(double* out, const double* xn) { int i = blockIdx.x * blockDim.x + threadIdx.x; double abar = 0.0, ye = 0.0; for (int n = 0; n < 4; ++n) { abar += xn[4 * i + n] / aion[n]; ye += xn[4 * i + n] * zion[n] / aion[n]; } out[i] = ye / abar; }
}
