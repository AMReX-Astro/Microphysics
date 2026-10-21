#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// burner/integrator pattern: a whole struct by non-const reference through two
// template layers, with an Array1D member updated via its accessor.
namespace probe {
template <class T, int XLO, int XHI> struct Array1D { T arr[XHI - XLO + 1]; __device__ T& operator()(int i) { return arr[i - XLO]; } };
struct Dvode { double t; int n_step; Array1D<double, 1, 14> y; };
template <class D> __device__ void advance(D& vs, double h) { for (int n = 1; n <= 14; ++n) vs.y(n) += h * n; vs.t += h; ++vs.n_step; }
template <class D> __device__ void integrate(D& vs) { for (int s = 0; s < 3; ++s) advance(vs, 0.5); }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; Dvode vs{}; vs.y(1) = in[i]; integrate(vs); out[i] = vs.y(1) + vs.t + vs.n_step; }
}
