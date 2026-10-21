#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// amrex::constexpr_for with a generic lambda (expected rejected).
namespace probe {
template <int I> struct IC { static constexpr int value = I; constexpr operator int() const { return I; } };
template <int B, int E, class F> __device__ void cfor(F&& f) { if constexpr (B < E) { f(IC<B>{}); cfor<B + 1, E>(f); } }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; double s = 0.0; cfor<1, 14>([&](auto n) { constexpr int m = n; s += m * in[i]; }); out[i] = s; }
}
