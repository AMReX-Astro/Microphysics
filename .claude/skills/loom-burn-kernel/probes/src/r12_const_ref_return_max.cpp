#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// amrex::max(const T&, const T&) -> const T& (reference return; expected rejected).
namespace probe {
template <class T> __device__ constexpr const T& mx(const T& a, const T& b) { return (a < b) ? b : a; }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; double a = in[i]; double b = 1.0; out[i] = mx(a, b); }
}
