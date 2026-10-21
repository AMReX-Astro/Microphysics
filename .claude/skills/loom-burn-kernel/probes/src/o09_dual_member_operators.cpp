#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// autodiff::dual style overloading with *member* operators (instance calls, hrx #2040).
namespace probe {
struct Dual {
  double v; double d;
  __device__ Dual operator*(Dual b) const { return Dual{v * b.v, d * b.v + v * b.d}; }
  __device__ Dual operator+(Dual b) const { return Dual{v + b.v, d + b.d}; }
  __device__ Dual& operator+=(Dual b) { v += b.v; d += b.d; return *this; }
};
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; Dual t{in[i], 1.0}; Dual r = t * t + t; r += t; out[i] = r.v + r.d; }
}
