#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// o08 by reference: auto& [a, b] = p must alias p's members.
namespace probe {
struct Pair { double first; double second; };
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; Pair p{in[i], 2.0}; auto& [a, b] = p; a = a * 3.0; b += a; out[i] = p.first + 10.0 * p.second; }
}
