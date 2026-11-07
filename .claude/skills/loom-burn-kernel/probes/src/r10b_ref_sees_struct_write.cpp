#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// r10 both ways: write through the struct, read through the ref; then write
// through the ref and copy the whole struct.
namespace probe {
struct S { double x; double y; };
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; S s{}; double& r = s.x; s.x = 2.0 * in[i]; double a = r; r = a + 5.0; S c = s; out[i] = a + 100.0 * c.x; }
}
