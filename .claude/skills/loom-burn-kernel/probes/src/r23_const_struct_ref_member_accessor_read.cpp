#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// Array1D member read through its const accessor via a const struct reference (integrator_to_burn(const dvode_t& vs, ...): vs.y(n)).
namespace probe {
template <class T, int XLO, int XHI> struct Array1D { T arr[XHI + 1]; __device__ T& operator()(int i) { return arr[i]; } __device__ const T& operator()(int i) const { return arr[i]; } };
struct Dvode { double t; int n; Array1D<double, 1, 14> y; };
struct Burn { double T; double xn[13]; };
__device__ static void to_burn(const Dvode& vs, Burn* s) { for (int n = 1; n <= 13; ++n) s->xn[n - 1] = vs.y(n); s->T = vs.t; }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; Dvode vs{}; for (int n = 1; n <= 14; ++n) vs.y(n) = in[i] * n; Burn s{}; to_burn(vs, &s); out[i] = s.xn[0] + s.T; }
}
