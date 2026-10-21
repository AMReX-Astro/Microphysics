#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// the referenced struct's own member functions write its Array1D member and scalar fields (`vs.add_y(n, v); vs.tick(h)`).
namespace probe {
template <class T, int XLO, int XHI> struct Array1D { T arr[XHI - XLO + 1]; __device__ T& operator()(int i) { return arr[i - XLO]; } };
struct Dvode { double t; int n_step; Array1D<double, 1, 14> y; __device__ void add_y(int n, double v) { y(n) += v; } __device__ void tick(double h) { t += h; ++n_step; } };
template <class D> __device__ void advance(D& vs, double h) { for (int n = 1; n <= 14; ++n) vs.add_y(n, h * n); vs.tick(h); }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; Dvode vs{}; vs.y(1) = in[i]; advance(vs, 0.5); out[i] = vs.y(1) + vs.t + vs.n_step; }
}
