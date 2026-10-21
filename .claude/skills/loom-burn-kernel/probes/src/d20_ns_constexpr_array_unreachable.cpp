#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// namespace-scope constexpr array that the root never reads (header tables the burn doesn't use).
namespace probe {
constexpr double unused_tab[3] = {4.0, 12.0, 16.0};
inline double host_only(int n) { return unused_tab[n]; }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; out[i] = in[i] * 2.0; }
}
