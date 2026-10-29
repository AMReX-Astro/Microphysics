#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// whole-struct assignment through a mutable struct reference (fill_additional_rates: rates1 = intermediate_rates(index)).
namespace probe {
template <class T, int XLO, int XHI> struct Array1D { T arr[XHI + 1]; __device__ T& operator()(int i) { return arr[i]; } __device__ const T& operator()(int i) const { return arr[i]; } };
template <class T, int XLO, int XHI, int YLO, int YHI> struct Array2D { T arr[(XHI + 1) * (YHI + 1)]; __device__ T& operator()(int i, int j) { return arr[i * (YHI + 1) + j]; } __device__ const T& operator()(int i, int j) const { return arr[i * (YHI + 1) + j]; } };
struct Burn { double rho; double T; double xn[13]; int n_rhs; };
struct Dvode { double t; int n; Array1D<double, 1, 14> y; Array2D<double, 1, 14, 1, 14> jac; };
template <class T> struct rate_t { T fr; T rr; };
struct Dual { double val; double grad; };
__device__ static void mk(rate_t<double>& r, double x) { r = rate_t<double>{x, 2.0 * x}; }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; rate_t<double> r{}; mk(r, in[i]); out[i] = r.fr + r.rr; }
}
