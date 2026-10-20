#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// amrex::Array1D<Real,1,N> (1-based operator() returning T& / const T&)
// passed by reference to an RHS-like function: actual_rhs(burn_t&, Array1D& ydot).
namespace probe {
template <class T, int XLO, int XHI> struct Array1D {
  T arr[XHI - XLO + 1];
  __device__ T& operator()(int i) { return arr[i - XLO]; }
  __device__ const T& operator()(int i) const { return arr[i - XLO]; }
};
struct Burn { double rho; double T; double xn[13]; };
template <class BurnT> __device__ void rhs(const BurnT& s, Array1D<double, 1, 14>& ydot) {
  for (int n = 1; n <= 13; ++n) ydot(n) = -s.rho * s.xn[n - 1];
  ydot(14) = s.T;
}
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; Burn s{}; s.rho = in[i]; s.T = 2.0; for (int n = 0; n < 13; ++n) s.xn[n] = 0.1 * n; Array1D<double, 1, 14> y{}; rhs(s, y); const auto& cy = y; out[i] = cy(1) + cy(14) + y(13); }
}
