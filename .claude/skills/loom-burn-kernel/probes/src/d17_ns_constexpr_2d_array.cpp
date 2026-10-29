#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// namespace-scope constexpr 2-D table (rate fit coefficients) with two runtime indices.
namespace probe {
constexpr double coef[3][4] = {{1.0, 2.0, 3.0, 4.0}, {5.0, 6.0, 7.0, 8.0}, {9.0, 10.0, 11.0, 12.0}};
KERNEL void k(double* out, const double* in, const int* sp) { int i = blockIdx.x * blockDim.x + threadIdx.x; int r = sp[i] % 3; double s = 0.0; for (int j = 0; j < 4; ++j) s = s * in[i] + coef[r][j]; out[i] = s; }
}
