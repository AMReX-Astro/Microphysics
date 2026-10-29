#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// amrex::Array2D<Real,1,N,1,N> accessor + pivoted LU in place via reference (dgefa shape).
namespace probe {
template <class T, int N> struct Array2D {
  T arr[N * N];
  __device__ T& operator()(int i, int j) { return arr[(j - 1) * N + (i - 1)]; }
  __device__ const T& operator()(int i, int j) const { return arr[(j - 1) * N + (i - 1)]; }
};
template <int N> __device__ void scale(Array2D<double, N>& a, int* ipvt) {
  for (int k = 1; k <= N; ++k) { ipvt[k - 1] = k; double t = -1.0 / a(k, k); for (int j = k + 1; j <= N; ++j) a(j, k) *= t; }
}
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; Array2D<double, 4> a{}; for (int r = 1; r <= 4; ++r) for (int c = 1; c <= 4; ++c) a(r, c) = in[i] + r + 2 * c; int p[4] = {}; scale(a, p); out[i] = a(4, 1) + a(2, 2) + p[3]; }
}
