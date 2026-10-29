#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// 40-deep nested ternary (generated rate tables / pynucastro-style selections; hrx #2079).
namespace probe {
__device__ double sel(double t, double x) { return t < 1.0 ? x * 1.0 : (t < 2.0 ? x * 2.0 : (t < 3.0 ? x * 3.0 : (t < 4.0 ? x * 4.0 : (t < 5.0 ? x * 5.0 : (t < 6.0 ? x * 6.0 : (t < 7.0 ? x * 7.0 : (t < 8.0 ? x * 8.0 : (t < 9.0 ? x * 9.0 : (t < 10.0 ? x * 10.0 : (t < 11.0 ? x * 11.0 : (t < 12.0 ? x * 12.0 : (t < 13.0 ? x * 13.0 : (t < 14.0 ? x * 14.0 : (t < 15.0 ? x * 15.0 : (t < 16.0 ? x * 16.0 : (t < 17.0 ? x * 17.0 : (t < 18.0 ? x * 18.0 : (t < 19.0 ? x * 19.0 : (t < 20.0 ? x * 20.0 : (t < 21.0 ? x * 21.0 : (t < 22.0 ? x * 22.0 : (t < 23.0 ? x * 23.0 : (t < 24.0 ? x * 24.0 : (t < 25.0 ? x * 25.0 : (t < 26.0 ? x * 26.0 : (t < 27.0 ? x * 27.0 : (t < 28.0 ? x * 28.0 : (t < 29.0 ? x * 29.0 : (t < 30.0 ? x * 30.0 : (t < 31.0 ? x * 31.0 : (t < 32.0 ? x * 32.0 : (t < 33.0 ? x * 33.0 : (t < 34.0 ? x * 34.0 : (t < 35.0 ? x * 35.0 : (t < 36.0 ? x * 36.0 : (t < 37.0 ? x * 37.0 : (t < 38.0 ? x * 38.0 : (t < 39.0 ? x * 39.0 : (t < 40.0 ? x * 40.0 : (x * 41.0)))))))))))))))))))))))))))))))))))))))); }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; out[i] = sel(in[i], in[i] + 1.0); }
}
