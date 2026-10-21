#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// namespace-scope constexpr array of aggregates (rate-parameter table), runtime index, field reads.
namespace probe {
struct Rp { double a; double b; int n; };
constexpr Rp rp[3] = {{1.0, 2.0, 1}, {3.0, 4.0, 2}, {5.0, 6.0, 3}};
KERNEL void k(double* out, const double* in, const int* sp) { int i = blockIdx.x * blockDim.x + threadIdx.x; const int r = sp[i] % 3; out[i] = rp[r].a * in[i] + rp[r].b * rp[r].n; }
}
