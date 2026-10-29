#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// replacement: template-struct recursion with static member functions and a pointer accumulator.
namespace probe {
template <int B, int E> struct Unroll { __device__ static void run(double* s, double x) { if constexpr (B < E) { *s += B * x; Unroll<B + 1, E>::run(s, x); } } };
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; double s = 0.0; Unroll<1, 14>::run(&s, in[i]); out[i] = s; }
}
