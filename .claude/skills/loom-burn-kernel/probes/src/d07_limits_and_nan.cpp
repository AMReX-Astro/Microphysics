#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// numeric_limits<double>::max/quiet_NaN replacements via builtins and bit casts.
namespace probe {
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; const double nan = __builtin_bit_cast(double, 0x7ff8000000000000ULL); const double big = __DBL_MAX__; out[i] = in[i] > 0.0 ? big : nan; }
}
