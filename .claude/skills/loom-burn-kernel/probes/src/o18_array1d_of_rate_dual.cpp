#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// Array1D<rate_t<Dual>, 1, N> intermediate-rate storage: element assignment and copies through template-struct recursion.
namespace probe {
template <class T, int XLO, int XHI> struct Array1D { T arr[XHI + 1]; __device__ T& operator()(int i) { return arr[i]; } __device__ const T& operator()(int i) const { return arr[i]; } };
template <class T, int XLO, int XHI, int YLO, int YHI> struct Array2D { T arr[(XHI + 1) * (YHI + 1)]; __device__ T& operator()(int i, int j) { return arr[i * (YHI + 1) + j]; } __device__ const T& operator()(int i, int j) const { return arr[i * (YHI + 1) + j]; } };
struct Burn { double rho; double T; double xn[13]; int n_rhs; };
struct Dvode { double t; int n; Array1D<double, 1, 14> y; Array2D<double, 1, 14, 1, 14> jac; };
template <class T> struct rate_t { T fr; T rr; };
struct Dual { double val; double grad; };
__device__ static rate_t<Dual> mkr(double x, int k) { return rate_t<Dual>{Dual{x * k, 1.0}, Dual{x + k, 2.0}}; }
template <int rate> struct Fill { static __device__ void run(double x, Array1D<rate_t<Dual>, 1, 5>& ir) { if constexpr (rate <= 5) { ir(rate) = mkr(x, rate); Fill<rate + 1>::run(x, ir); } } };
template <int rate> struct Use { static __device__ void run(const Array1D<rate_t<Dual>, 1, 5>& ir, double* acc) { if constexpr (rate <= 5) { const rate_t<Dual> r = ir(rate); acc[0] += r.fr.val * r.rr.grad + r.fr.grad; Use<rate + 1>::run(ir, acc); } } };
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; Array1D<rate_t<Dual>, 1, 5> ir{}; Fill<1>::run(in[i], ir); double acc[1] = {0.0}; Use<1>::run(ir, acc); out[i] = acc[0]; }
}
