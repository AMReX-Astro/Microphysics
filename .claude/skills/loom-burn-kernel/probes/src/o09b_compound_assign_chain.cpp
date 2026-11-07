#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// o09 with operator+= results used as references: (r += t) += t, and a ref to the result.
namespace probe {
struct Dual {
  double v; double d;
  __device__ Dual operator*(Dual b) const { return Dual{v * b.v, d * b.v + v * b.d}; }
  __device__ Dual& operator+=(Dual b) { v += b.v; d += b.d; return *this; }
};
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; Dual t{in[i], 1.0}; Dual r = t * t; (r += t) += t; Dual& q = (r += r); q.d += 0.5; out[i] = r.v + 1000.0 * r.d; }
}
