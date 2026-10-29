#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// auto [fwd, rev] = rhs_term(...) returning a 2-field struct (std::pair replaced).
namespace probe {
struct Pair { double first; double second; };
__device__ Pair term(double x) { return Pair{x * 2.0, -x}; }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; auto [f, r] = term(in[i]); out[i] = f + r; }
}
