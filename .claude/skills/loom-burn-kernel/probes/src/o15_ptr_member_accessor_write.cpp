#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// same writes through a pointer (`vs->y(n) += ...; vs->t += h`), the pointer taken from a reference parameter.
namespace probe {
template <class T, int XLO, int XHI> struct Array1D { T arr[XHI - XLO + 1]; __device__ T& operator()(int i) { return arr[i - XLO]; } };
struct Dvode { double t; int n_step; Array1D<double, 1, 14> y; __device__ void add_y(int n, double v) { y(n) += v; } __device__ void tick(double h) { t += h; ++n_step; } };
template <class D> __device__ void advance_p(D* vs, double h) { for (int n = 1; n <= 14; ++n) vs->y(n) += h * n; vs->t += h; ++vs->n_step; }
template <class D> __device__ void advance(D& vs, double h) { advance_p(&vs, h); }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; Dvode vs{}; vs.y(1) = in[i]; advance(vs, 0.5); out[i] = vs.y(1) + vs.t + vs.n_step; }
}
