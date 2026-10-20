#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// const BurnParams& bound to *params, a kernel-argument pointer (burner(&s, dt, *params, tab)).
namespace probe {
template <class T, int XLO, int XHI> struct Array1D { T arr[XHI + 1]; __device__ T& operator()(int i) { return arr[i]; } __device__ const T& operator()(int i) const { return arr[i]; } };
template <class T, int XLO, int XHI, int YLO, int YHI> struct Array2D { T arr[(XHI + 1) * (YHI + 1)]; __device__ T& operator()(int i, int j) { return arr[i * (YHI + 1) + j]; } __device__ const T& operator()(int i, int j) const { return arr[i * (YHI + 1) + j]; } };
struct Burn { double rho; double T; double xn[13]; int n_rhs; };
struct Dvode { double t; int n; Array1D<double, 1, 14> y; Array2D<double, 1, 14, 1, 14> jac; };
template <class T> struct rate_t { T fr; T rr; };
struct Dual { double val; double grad; };
struct EP { double mintemp; int flag; };
struct IP { double rtol; double atol; int n; };
struct BP { EP eos; IP integ; };
__device__ static double use(const IP& ip, const EP& ep, double x) { return ep.flag ? x * ip.rtol + ip.atol : x - ep.mintemp; }
__device__ static double burner(const BP& p, double x) { const EP ep = p.eos; const IP ip = p.integ; return use(ip, ep, x); }
KERNEL void k(const BP* params, double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; out[i] = burner(*params, in[i]); }
}
