#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// dual_array<1,3> (sneut5, screening) as POD with a 3-gradient, free functions.
namespace probe {
struct Dual3 { double v; double g[3]; };
__device__ Dual3 dmul(Dual3 a, Dual3 b) { Dual3 r{}; r.v = a.v * b.v; for (int n = 0; n < 3; ++n) r.g[n] = a.g[n] * b.v + a.v * b.g[n]; return r; }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; Dual3 t{}; t.v = in[i]; t.g[0] = 1.0; Dual3 r = dmul(t, t); out[i] = r.v + r.g[0] + r.g[2]; }
}
