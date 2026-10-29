#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// autodiff::dual style operator overloading (expected rejected).
namespace probe {
struct Dual { double v; double d; };
__device__ Dual operator*(Dual a, Dual b) { return Dual{a.v * b.v, a.d * b.v + a.v * b.d}; }
__device__ Dual operator+(Dual a, Dual b) { return Dual{a.v + b.v, a.d + b.d}; }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; Dual t{in[i], 1.0}; Dual r = t * t + t; out[i] = r.v + r.d; }
}
