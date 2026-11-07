#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// l01 with two by-reference captures, a by-value capture, and the lambda called twice.
namespace probe {
template <int I> struct IC { static constexpr int value = I; constexpr operator int() const { return I; } };
template <int B, int E, class F> __device__ void cfor(F&& f) { if constexpr (B < E) { f(IC<B>{}); cfor<B + 1, E>(f); } }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; double s = 0.0, p = 1.0; const double x = in[i];
  auto body = [&, x](auto n) { constexpr int m = n; s += m * x; p *= (1.0 + 0.125 * m); };
  cfor<1, 6>(body); cfor<2, 4>(body); out[i] = s + p; }
}
