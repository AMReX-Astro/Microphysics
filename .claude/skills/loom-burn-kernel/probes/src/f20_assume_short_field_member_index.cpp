#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// `const int La = vs->L` (a short field), assumed in [1, 6], indexes an Array1D member through a pointer (dvset el(La)).
namespace probe {
template <class T, int XLO, int XHI> struct Array1D { T arr[XHI + 1]; __device__ T& operator()(int i) { return arr[i]; } __device__ const T& operator()(int i) const { return arr[i]; } };
struct Dvode { short L; short NQ; Array1D<double, 1, 6> el; Array1D<double, 1, 13> tau; };
__device__ static void set(Dvode* vs) { const int La = vs->L; __builtin_assume(1 <= La && La <= 6); vs->el(1) = vs->el(La) * 2.0; }
KERNEL void k(double* out, const double* in, const short* l_in) { int i = blockIdx.x * blockDim.x + threadIdx.x; Dvode vs{}; for (int n = 1; n <= 6; ++n) vs.el(n) = in[i] + n; vs.L = l_in[i]; set(&vs); out[i] = vs.el(1); }
}
