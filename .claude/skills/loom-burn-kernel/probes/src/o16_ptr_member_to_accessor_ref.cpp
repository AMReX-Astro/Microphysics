#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// Array1D member reached through a pointer, passed to a mutable Array1D& and written (rhs(..., vs->savf)).
namespace probe {
template <class T, int XLO, int XHI> struct Array1D { T arr[XHI + 1]; __device__ T& operator()(int i) { return arr[i]; } __device__ const T& operator()(int i) const { return arr[i]; } };
template <class T, int XLO, int XHI, int YLO, int YHI> struct Array2D { T arr[(XHI + 1) * (YHI + 1)]; __device__ T& operator()(int i, int j) { return arr[i * (YHI + 1) + j]; } __device__ const T& operator()(int i, int j) const { return arr[i * (YHI + 1) + j]; } };
struct Burn { double rho; double T; double xn[13]; int n_rhs; };
struct Dvode { double t; int n; Array1D<double, 1, 14> y; Array2D<double, 1, 14, 1, 14> jac; };
template <class T> struct rate_t { T fr; T rr; };
struct Dual { double val; double grad; };
__device__ static void upd(Array1D<double, 1, 14>& y, double h) { for (int n = 1; n <= 14; ++n) y(n) += h * n; }
__device__ static void adv(Dvode* vs, double h) { upd(vs->y, h); vs->t += h; }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; Dvode vs{}; vs.y(1) = in[i]; adv(&vs, 0.5); out[i] = vs.y(1) + vs.y(14) + vs.t; }
}
