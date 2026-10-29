#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// field of a const struct reference passed by value to a call that also takes a double& out-parameter
// (`f(s.eta, o)` with `const S& s` and an initialized local o; aprox19 ecapnuc(burn_state.eta, ...)).
namespace probe {
struct S { double eta; double T; };
static void set(const double x, double& o) { o = 2.0 * x; }
static double g(const S& s) { double o = 0.0; set(s.eta, o); return o; }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; S s{}; s.eta = in[i]; out[i] = g(s); }
}
