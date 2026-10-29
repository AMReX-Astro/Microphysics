#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// non-operator member function on a local object (VODE / Array helpers).
namespace probe {
struct Acc { double s; __device__ void add(double x) { s += x; } };
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; Acc a{}; a.add(in[i]); out[i] = a.s; }
}
