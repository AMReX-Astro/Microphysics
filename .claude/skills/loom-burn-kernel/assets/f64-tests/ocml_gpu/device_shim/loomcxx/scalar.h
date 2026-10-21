// Device (hipcc, AMDGPU) stand-in for <loomcxx/scalar.h>, for compiling the
// f64 OCML recipes (../../../../kernel/f64_math.H, f64_math_extra.H) as HIP
// __device__ functions. See ../../README.md.
//
// Each loom::scalar operation the recipes use is mapped to the AMDGPU
// operation Loom is expected to emit for it:
//
//   fmaf(double)              __builtin_fma        llvm.fma.f64   v_fma_f64
//   truncf(double)            __builtin_trunc      llvm.trunc.f64 v_trunc_f64
//   roundevenf(double/float)  __builtin_roundeven  llvm.roundeven v_rndne_f64/f32
//   ldexpf(double)            __builtin_ldexp      llvm.ldexp     v_ldexp_f64
//   frexpf(double)            __builtin_frexp      llvm.frexp     v_frexp_mant_f64,
//                                                                 v_frexp_exp_i32_f64
//   approximate::log2f(float) __builtin_amdgcn_logf  llvm.amdgcn.log.f32  v_log_f32
//   approximate::exp2f(float) __builtin_amdgcn_exp2f llvm.amdgcn.exp2.f32 v_exp_f32
//
// fma, trunc, roundeven and ldexp are exact IEEE operations and frexp is an
// exact decomposition; any correct lowering gives the same bits. These are
// also the builtins OCML itself calls for BUILTIN_FMA/RINT/FLDEXP/FREXP_*.
// The approximate f32 log2/exp2 (only f64_math.H's cbrt uses them, for its
// seed) are the uncertain part: OCML's cbrtD.cl uses
// BUILTIN_AMDGPU_LOG2_F32/BUILTIN_AMDGPU_EXP2_F32, i.e. the raw instructions,
// which is what the default variant here emits. If Loom instead emits
// llvm.log2.f32/llvm.exp2.f32 with the afn flag, the AMDGPU backend still
// selects v_log_f32/v_exp_f32 but wraps them in a denormal-input range
// check; the cbrt arguments are in [0.5, 2) (log2) and about [-0.34, 0.34]
// (exp2), where that check is a no-op, so the bits are the same. Select a
// different lowering with -DOCMLGPU_APPROX_VARIANT=N to see whether it
// matters:
//   0 (default)  __builtin_amdgcn_logf / __builtin_amdgcn_exp2f
//   1            __ocml_native_log2_f32 / __ocml_native_exp2_f32 (device-libs)
//   2            __builtin_log2f / __builtin_exp2f without afn (the accurate
//                llvm.log2.f32 / llvm.exp2.f32 lowering; NOT what OCML uses,
//                for diagnosis only)
//
// The other operations the host shim provides (sqrtf, floorf) are not used by
// these recipes; they are mapped to the
// corresponding llvm intrinsics for completeness (llvm.sqrt.f64 is the
// correctly rounded lowering, v_sqrt_f64 plus refinement).
//
// Everything else in the recipes (f64 add/sub/mul, IEEE f64 division, which
// AMDGPU lowers to the correctly rounded v_div_scale/v_rcp/v_fma/v_div_fmas/
// v_div_fixup sequence, compares, selects, f64<->f32 and int conversions,
// integer and bit operations) is plain C++ that hipcc lowers directly.
#pragma once

#ifndef OCMLGPU_APPROX_VARIANT
#define OCMLGPU_APPROX_VARIANT 0
#endif

#if OCMLGPU_APPROX_VARIANT == 1
extern "C" __device__ float __ocml_native_log2_f32(float);
extern "C" __device__ float __ocml_native_exp2_f32(float);
#endif

namespace loom::scalar {
__device__ inline double fmaf(double a, double b, double c) { return __builtin_fma(a, b, c); }
__device__ inline float fmaf(float a, float b, float c) { return __builtin_fmaf(a, b, c); }
__device__ inline double truncf(double x) { return __builtin_trunc(x); }
__device__ inline float truncf(float x) { return __builtin_truncf(x); }
__device__ inline double absf(double x) { return __builtin_fabs(x); }
__device__ inline double copysignf(double x, double y) { return __builtin_copysign(x, y); }
__device__ inline double sqrtf(double x) { return __builtin_sqrt(x); }
__device__ inline double floorf(double x) { return __builtin_floor(x); }
__device__ inline double roundevenf(double x) { return __builtin_roundeven(x); }
__device__ inline float roundevenf(float x) { return __builtin_roundevenf(x); }
__device__ inline float absf(float x) { return __builtin_fabsf(x); }
__device__ inline double ldexpf(double x, int n) { return __builtin_ldexp(x, n); }
template <class Float> struct frexpf_result {
    Float fraction;
    int exponent;
};
__device__ inline frexpf_result<double> frexpf(double x) {
    int e;
    const double m = __builtin_frexp(x, &e);
    return {m, e};
}

namespace approximate {
__device__ inline float log2f(float x) {
#if OCMLGPU_APPROX_VARIANT == 1
    return __ocml_native_log2_f32(x);
#elif OCMLGPU_APPROX_VARIANT == 2
    return __builtin_log2f(x);
#else
    return __builtin_amdgcn_logf(x);
#endif
}
__device__ inline float exp2f(float x) {
#if OCMLGPU_APPROX_VARIANT == 1
    return __ocml_native_exp2_f32(x);
#elif OCMLGPU_APPROX_VARIANT == 2
    return __builtin_exp2f(x);
#else
    return __builtin_amdgcn_exp2f(x);
#endif
}
} // namespace approximate
} // namespace loom::scalar
