#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// eos.H reset_inputs: switch on a runtime eos_input_t (expected rejected).
namespace probe {
enum class In : int { rt = 1, rh = 2, re = 5 };
KERNEL void k(double* out, const double* in, const int* mode) { int i = blockIdx.x * blockDim.x + threadIdx.x; double r = 0.0;
  switch (static_cast<In>(mode[i])) { case In::rt: r = in[i]; break; case In::re: r = 2.0 * in[i]; break; default: r = -1.0; break; } out[i] = r; }
}
