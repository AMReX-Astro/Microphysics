#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// uninitialized scalar local bound to a double& out-parameter (`double o; set(x, o);`, aprox19 ecapnuc's tmp_fr).
namespace probe {
static void set(const double x, double& o) { o = 2.0 * x; }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; double o; set(in[i], o); out[i] = o; }
}
