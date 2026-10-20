#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// dgefa on an Array2D member through a pointer, runtime bounds with __builtin_assume (rhs_dgefa(vs->jac)).
namespace probe {
template <class T, int XLO, int XHI> struct Array1D { T arr[XHI + 1]; __device__ T& operator()(int i) { return arr[i]; } __device__ const T& operator()(int i) const { return arr[i]; } };
template <class T, int XLO, int XHI, int YLO, int YHI> struct Array2D { T arr[(XHI + 1) * (YHI + 1)]; __device__ T& operator()(int i, int j) { return arr[i * (YHI + 1) + j]; } __device__ const T& operator()(int i, int j) const { return arr[i * (YHI + 1) + j]; } };
struct Burn { double rho; double T; double xn[13]; int n_rhs; };
struct Dvode { double t; int n; Array1D<double, 1, 14> y; Array2D<double, 1, 14, 1, 14> jac; };
template <class T> struct rate_t { T fr; T rr; };
struct Dual { double val; double grad; };
__device__ static int dgefa(Array2D<double, 1, 14, 1, 14>& a, int nn) { int info = 0; for (int k = 1; k <= nn - 1; ++k) { __builtin_assume(1 <= k && k <= 13); if (a(k, k) == 0.0) { info = k; } else { double t = -1.0 / a(k, k); for (int j = k + 1; j <= nn; ++j) { __builtin_assume(2 <= j && j <= 14); a(j, k) *= t; } for (int j = k + 1; j <= nn; ++j) { __builtin_assume(2 <= j && j <= 14); t = a(k, j); for (int ii = k + 1; ii <= nn; ++ii) { __builtin_assume(2 <= ii && ii <= 14); a(ii, j) += t * a(ii, k); } } } } return info; }
__device__ static void adv(Dvode* vs, int nn) { vs->n = dgefa(vs->jac, nn); }
KERNEL void k(double* out, const double* in, int nn) { int i = blockIdx.x * blockDim.x + threadIdx.x; Dvode vs{}; for (int r = 1; r <= 14; ++r) for (int c = 1; c <= 14; ++c) vs.jac(r, c) = in[i] + r + (r == c ? 20.0 : 0.5 * c); adv(&vs, nn); out[i] = vs.jac(14, 14) + vs.n; }
}
