#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// const Burn& bound to *p, p a mutable Burn* (burner -> network: actual_rhs(*state, ...)).
namespace probe {
template <class T, int XLO, int XHI> struct Array1D { T arr[XHI + 1]; __device__ T& operator()(int i) { return arr[i]; } __device__ const T& operator()(int i) const { return arr[i]; } };
template <class T, int XLO, int XHI, int YLO, int YHI> struct Array2D { T arr[(XHI + 1) * (YHI + 1)]; __device__ T& operator()(int i, int j) { return arr[i * (YHI + 1) + j]; } __device__ const T& operator()(int i, int j) const { return arr[i * (YHI + 1) + j]; } };
struct Burn { double rho; double T; double xn[13]; int n_rhs; };
struct Dvode { double t; int n; Array1D<double, 1, 14> y; Array2D<double, 1, 14, 1, 14> jac; };
template <class T> struct rate_t { T fr; T rr; };
struct Dual { double val; double grad; };
__device__ static double rd(const Burn& s) { return s.rho * s.T + s.xn[3]; }
__device__ static void step(Burn* p) { p->T = rd(*p); ++p->n_rhs; }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; Burn s{}; s.rho = in[i]; s.T = 2.0; step(&s); out[i] = s.T + s.n_rhs; }
}
