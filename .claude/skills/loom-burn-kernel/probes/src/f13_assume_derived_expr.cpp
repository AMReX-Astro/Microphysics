#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// assumption over a derived expression (`__builtin_assume(j + 1 <= n)`) guarding a private-array index.
namespace probe {
KERNEL void k(double* out, const double* in, const int* n_in) { int i = blockIdx.x * blockDim.x + threadIdx.x; double a[15] = {}; for (int n = 0; n < 15; ++n) a[n] = in[i] + n; int n = n_in[i]; __builtin_assume(1 <= n && n <= 14); double s = 0.0; for (int j = 1; j < n; ++j) { __builtin_assume(1 <= j && j + 1 <= n); s += a[j + 1] - a[j]; } out[i] = s; }
}
