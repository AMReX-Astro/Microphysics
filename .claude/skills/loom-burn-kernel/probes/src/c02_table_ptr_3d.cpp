#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// helmholtz f[jmax][imax][9] flattened, indexed from a log10-free hash; clamp like amrex::Clamp.
namespace probe {
constexpr int imax = 541; constexpr int jmax = 201;
__device__ int clampi(int v, int lo, int hi) { return v < lo ? lo : (hi < v ? hi : v); }
KERNEL void k(const double* f, double* out, const double* in) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int jat = clampi(static_cast<int>(in[2*i] * 50.0) + 1, 1, jmax - 1) - 1;
  int iat = clampi(static_cast<int>(in[2*i+1] * 30.0) + 1, 1, imax - 1) - 1;
  double fi[36] = {};
  for (int n = 0; n < 9; ++n) { fi[n] = f[(static_cast<long>(jat) * imax + iat) * 9 + n]; fi[9+n] = f[(static_cast<long>(jat+1) * imax + iat) * 9 + n]; }
  double s = 0.0; for (int n = 0; n < 18; ++n) s += fi[n];
  out[i] = s;
}
}
