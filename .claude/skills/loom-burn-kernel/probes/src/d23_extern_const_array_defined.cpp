#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// extern const array declared in a header and defined in the same TU (strong definition, #2163).
namespace probe {
extern const double aion[3];
extern const double aion[3] = {4.0, 12.0, 16.0};
KERNEL void k(double* out, const double* in, const int* sp) { int i = blockIdx.x * blockDim.x + threadIdx.x; out[i] = in[i] / aion[sp[i] % 3]; }
}
