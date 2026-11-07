// Usage: gpu_driver <code object> <kernel symbol> <nargs 2|3> <out file>
#include <hip/hip_runtime.h>
#include <cstdlib>
#include "common.H"
#define CK(x) do { hipError_t e_ = (x); if (e_ != hipSuccess) { std::fprintf(stderr, "%s: %s\n", #x, hipGetErrorString(e_)); std::exit(2); } } while (0)
int main(int argc, char** argv) {
    if (argc != 5) { std::fprintf(stderr, "usage\n"); return 2; }
    const int nargs = std::atoi(argv[3]);
    double in[N], out[N]; int sp[N];
    fill_inputs(in, sp);
    for (int i = 0; i < N; ++i) out[i] = -1.0;
    double *in_d, *out_d; int* sp_d;
    CK(hipMalloc(&in_d, sizeof in)); CK(hipMalloc(&out_d, sizeof out)); CK(hipMalloc(&sp_d, sizeof sp));
    CK(hipMemcpy(in_d, in, sizeof in, hipMemcpyHostToDevice));
    CK(hipMemcpy(out_d, out, sizeof out, hipMemcpyHostToDevice));
    CK(hipMemcpy(sp_d, sp, sizeof sp, hipMemcpyHostToDevice));
    hipModule_t mod; hipFunction_t fn;
    CK(hipModuleLoad(&mod, argv[1]));
    CK(hipModuleGetFunction(&fn, mod, argv[2]));
    double* a_out = out_d; const double* a_in = in_d; const int* a_sp = sp_d;
    void* args3[] = {&a_out, &a_in, &a_sp};
    void* args2[] = {&a_out, &a_in};
    CK(hipModuleLaunchKernel(fn, 1, 1, 1, N, 1, 1, 0, nullptr, nargs == 3 ? args3 : args2, nullptr));
    CK(hipDeviceSynchronize());
    CK(hipMemcpy(out, out_d, sizeof out, hipMemcpyDeviceToHost));
    dump(argv[4], out);
}
