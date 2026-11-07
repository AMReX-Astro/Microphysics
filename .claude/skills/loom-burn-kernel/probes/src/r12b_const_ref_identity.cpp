#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// r12 identity: keep the returned const T&, then change the selected variable.
// A by-value fallback returns the old value.
namespace probe {
template <class T> __device__ constexpr const T& mx(const T& a, const T& b) { return (a < b) ? b : a; }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; double a = in[i]; double b = 1.0; const double& m = mx(a, b); a = 7.0 * a; b = -3.0; out[i] = m; }
}
