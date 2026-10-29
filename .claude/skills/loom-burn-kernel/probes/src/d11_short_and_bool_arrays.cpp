#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// Array1D<short> pivot and bool flags in a struct, int<->short conversions.
namespace probe {
struct P { short piv[14]; bool used[14]; };
KERNEL void k(double* out, const double* in, const int* idx) { int i = blockIdx.x * blockDim.x + threadIdx.x; P p{}; for (int n = 0; n < 14; ++n) { p.piv[n] = static_cast<short>(13 - n); p.used[n] = (n % 2) == 0; } int j = p.piv[idx[i] % 14]; out[i] = in[i] + j + (p.used[j] ? 1.0 : 0.0); }
}
