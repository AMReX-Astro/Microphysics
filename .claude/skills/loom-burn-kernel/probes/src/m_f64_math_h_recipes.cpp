#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// OCML f64 transcriptions: the skill's assets/kernel/f64_math.H, staged as f64_math.h by
// probe.sh (F64_INC). exp, log, cbrt only (pow/log10/tanh are in f64_math_extra.H).
// f64_math.h expects DEVICE and lives in namespace chemistry::f64_math.
#define DEVICE __device__
#include "f64_math.h"
namespace probe {
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x;
  out[i] = chemistry::f64_math::exp(in[i]) + chemistry::f64_math::log(in[i] + 2.0) + chemistry::f64_math::cbrt(in[i]); }
}
