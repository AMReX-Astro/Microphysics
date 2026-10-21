# GPU bit-for-bit check of the OCML f64 recipes

`../../kernel/f64_math.H` (`exp`, `log`, `cbrt`) and
`../../kernel/f64_math_extra.H` (`pow`, `log10`, `tanh`) transcribe AMD OCML
(llvm-project `amd/device-libs` d366fa84f3fd) because Loom cannot lower these
f64 functions on AMDGPU. A Loom burn must reproduce the stock HIP build bit
for bit, and the stock build calls OCML through `std::pow/log10/tanh/exp/log/
cbrt` in device code. The host test in `..` checks the recipes against MPFR;
this test runs on the GPU and compares the recipes' bits with hipcc's own
results for the same inputs. It is network-neutral: the sweeps cover the
argument ranges of Microphysics' rate, EOS and VODE code in general.

## What it checks

For every input, three results are computed on the GPU:

- `recipe`: `chemistry::f64_math::<fn>` compiled as a `__device__` function,
  with `device_shim/loomcxx/scalar.h` mapping each `loom::scalar` op to the
  AMDGPU instruction Loom is expected to emit (`fmaf` -> `v_fma_f64`,
  `truncf` -> `v_trunc_f64`, `approximate::log2f/exp2f` -> `v_log_f32`/
  `v_exp_f32`);
- `stock`: `std::<fn>` in device code, as the stock Microphysics HIP build
  writes it;
- `direct`: `__ocml_<fn>_f64` by name (informational).

The host compares `recipe` vs `stock` bit patterns per sweep: special and
edge values (signed zeros, subnormals, infinities, NaNs with payloads,
overflow/underflow thresholds), random bit patterns, and function-specific
ranges (pow with run-time and compile-time-constant exponents, base 10,
VODE's `pow(x, 1/L)` and `pow(x, L)`, near-1 bases, integer and half-integer
exponents; log/log10/cbrt/exp/tanh over their domains). Same generators and
seed as `../test_ocml.cpp`.

Policy: results where both sides are NaN with different payloads are
counted separately (`nan_payload_diff`) and do not count as `differ`; NaN
payloads are not part of the bit-for-bit requirement. Extra diagnostic
lines: `OCMLGPU pow(x,1)==x` (OCML pow has no `y == 1` case) and
`OCMLGPU-DIAG rcp/...` (IEEE division vs an emulated `MATH_FAST_RCP`).

## Build and run (on the GPU host)

Ship the skill's `assets/` directory (it needs `kernel/` and
`f64-tests/ocml_gpu/` at the same relative positions), then on the host:

```bash
cd assets/f64-tests/ocml_gpu
./build_and_run.sh          # 1048576 samples per random sweep; hipcc from PATH or $ROCM_PATH/bin (default /opt/rocm)
./build_and_run.sh 4194304                                     # more samples
OFFLOAD_ARCH=gfx942 ./build_and_run.sh                         # another target
OCMLGPU_APPROX_VARIANT=1 ./build_and_run.sh                    # other f32 log2/exp2 lowering
./build_and_run.sh -h                                          # all options
```

It compiles with `hipcc -std=c++20 -O3 -ffp-contract=off
--offload-arch=$OFFLOAD_ARCH` (the stock reference flags), writes to
`$BUILD_DIR` (default `${TMPDIR:-/tmp}/loom_ocml_gpu`), dumps the device
assembly with instruction counts, and runs under `timeout` (`RUN_TIMEOUT`,
default 300 s). Last line: `PASS` only if every sweep is bit-identical
including NaN payloads; otherwise `FAIL: only NaN payloads differ` (acceptable
under the policy above) or `FAIL: some sweep has non-NaN bit differences`.
Exit status 0 on PASS, 1 on FAIL, 2 on a HIP error.

## Recorded result (gfx90a)

Bit-identical to hipcc's OCML for `pow`, `log10`, `tanh`, `exp`, `log` and
`cbrt` over about 1M samples per sweep. The only differences were 3 `cbrt`
results with different NaN payloads (`nan_payload_diff=3`), which the policy
ignores, so the run ends in `FAIL: only NaN payloads differ`.

2026-10-21, after the recipes switched to native rint/ldexp/frexp (OCML's
own builtins; only `MATH_FAST_RCP`/`MATH_FAST_DIV` remain IEEE division):
the same result at 4,194,304 samples per sweep (bit-identical everywhere,
the same 3 `cbrt` NaN payloads), and the `rcp` diagnostic found no input
where IEEE `1/x` differs from `v_rcp_f64` + 2 Newton steps.
