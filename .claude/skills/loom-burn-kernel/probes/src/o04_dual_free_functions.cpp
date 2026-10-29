#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// Proposed rewrite: dual as POD + free functions, by-value args and returns, same op order.
namespace probe {
struct Dual { double v; double d; };
__device__ Dual dmul(Dual a, Dual b) { const double v = a.v * b.v; const double d = a.d * b.v + a.v * b.d; return Dual{v, d}; }
__device__ Dual dadd(Dual a, Dual b) { return Dual{a.v + b.v, a.d + b.d}; }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; Dual t{in[i], 1.0}; Dual r = dadd(dmul(t, t), t); out[i] = r.v + r.d; }
}
