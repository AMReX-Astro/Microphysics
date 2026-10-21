# Host validation of the OCML f64 pow, log10 and tanh recipes

`../kernel/f64_math_extra.H` transcribes AMD OCML's binary64 `pow`, `log10` and
`tanh` for the Loom kernel, which cannot lower f64 pow/log10/tanh on AMDGPU.
The goal is that the Loom kernel computes the same bits as the HIP reference
build, where hipcc turns `std::pow(double, double)`, `std::log10(double)` and
`std::tanh(double)` into `__ocml_pow_f64`, `__ocml_log10_f64` and
`__ocml_tanh_f64`.

## OCML source revision

llvm-project `amd/device-libs` at commit **d366fa84f3fd** (the same revision
as `f64_math.H`), files in `amd/device-libs/ocml/src/`:

| Function | Files |
|---|---|
| `pow` | `powD.cl`, `powD_base.h` (`COMPILING_POW`), `eplnD.cl`, `expepD.cl`, `expD_base.h` (through `f64_math.H`'s `exp`), `ep.h` |
| `log10` | `log10D.cl`, `logD_base.h` (`COMPILING_LOG10`), `ep.h` |
| `tanh` | `tanhD.cl`, `epexpepD.cl`, `ep.h` |

Build controls that select the transcribed paths:

- `EXTRA_ACCURACY` is not defined anywhere in device-libs, so `expep` is
  `exp(x.hi)` followed by `fma(z, x.lo, z)`.
- HIP without `-ffast-math`/`-ffinite-math-only` links
  `__oclc_finite_only_opt = 0` and `__oclc_unsafe_math_opt = 0`, so
  `log10` keeps its non-finite branch. None of these functions read
  `DAZ_OPT()` (an f32 control); gfx90a keeps f64 denormals.
- OCML is compiled as OpenCL C (`-cl-std=CL2.0`), where contraction within a
  source expression is on. None of the transcribed expressions has a
  multiply feeding an add in the same expression (every fused operation is an
  explicit `BUILTIN_FMA_F64`/`MATH_MAD`), so contraction cannot change them.
  The kernel's own `-ffp-contract=off` does not reach OCML.

`pown` is not needed: nothing on the first port's (aprox13 + helmholtz + VODE) path calls
`amrex::Math::pown`, and the only integer-exponent `std::pow` call
(`std::pow(vstate.H / vstate.tau(2), vstate.L)` in `vode_dvstep.H`, `L` is
`short`) resolves on HIP to the `__HIP_OVERLOAD2` template (exact match beats
the promotion to `pow(double, int)`), which calls
`pow((double)x, (double)L)`, i.e. `__ocml_pow_f64`. Translate it as
`mpl::math::pow(x, (Real)L)`.

## Hardware-dependent steps

- **`MATH_FAST_RCP`** (`ep.h` `div`, `rcp`, `fdiv` in `epln`, `log10`, `tanh`):
  OCML refines the hardware `v_rcp_f64` with two fma Newton steps. The
  recipes use IEEE division, as `f64_math.H` does: Loom exposes no raw
  `v_rcp_f64` (probes m_f64_div_afn, m_f64_divide_reciprocal,
  m_f64_amdgcn_rcp). This is the only substitution left; rint, ldexp and
  frexp use the same native instructions as OCML (hrx-system #2131). Two Newton steps from
  `v_rcp_f64` give the correctly rounded reciprocal in essentially every case;
  the GPU bit check (`ocml_gpu/`) found no difference on gfx90a.
- No f32 approximate instructions (`v_log_f32`, `v_exp_f32`) are used by
  these three functions. (`f64_math.H`'s `cbrt` uses them; the host shim
  models them with correctly rounded host functions, so a host test of
  `cbrt` would be approximate for that step.)
- NaN payloads: OCML returns a propagated or canonical NaN; the recipes
  return the quieted input (log10, tanh, keeping its sign) or the canonical
  quiet NaN (pow). LLVM does not guarantee NaN payloads either. NaN sign and
  payload bits don't matter in any comparison (user policy): any NaN matches
  any NaN.

All other operations used (f64 add/sub/mul/div, fma, trunc, integer and bit
operations) are exact IEEE operations, so the host build computes the same
bits Loom will.

## Running

```bash
bash <skill>/assets/f64-tests/run.sh 2000000
# samples per random sweep; default 1000000. Recipes are taken from ../kernel/
```

`run.sh` compiles `test_ocml.cpp` with `clang++ -std=c++20 -ffp-contract=off`
at `-O2` and `-O0` against `host_shim/loomcxx/scalar.h` (host stand-ins for
the `loom::scalar` operations) and MPFR (macOS: `brew install mpfr`; Ubuntu:
`apt install libmpfr-dev libgmp-dev clang`; override the prefixes with
`MPFR_PREFIX`/`GMP_PREFIX` and the compiler with `CXX`; without a `clang++`
on PATH it falls back to ROCm's `$ROCM_PATH/llvm/bin/clang++`), puts binaries in `$BUILD_DIR` (default
`$TMPDIR/loom_ocml_test`), runs the sweeps (pow, log10, tanh from
`f64_math_extra.H`, and since 2026-10-19 also exp, log, cbrt from
`f64_math.H` directly), and always checks that the `-O0` and `-O2` builds
produce the same output bits (hash), even after a failed sweep. It exits
non-zero on a special-case mismatch, any error above 2 ulp, or differing
hashes. The hash depends on the C++ standard library's random distributions
(libc++ vs libstdc++), so compare hashes only between builds on one machine.

Each sweep reports the recipe's ulp distance from the correctly rounded
result (MPFR, 53 bits, IEEE subnormals) as a histogram, the macOS libm's
distance from the same reference, and the recipe-vs-libm distance. Results
that are NaN, infinite or zero must equal libm's bits (any NaN matches any
NaN). The pow special grid (24 x 29 inputs: signed zeros, smallest
subnormal, +-0.5, +-1, +-2, +-3, +-2.5, 10, 1 +- ulp, +-DBL_MAX, +-inf,
+-NaN against signed zeros, +-0.5, +-1, +-2, +-3, +-2.5, +-2^53, odd
integers near 2^53, 1023.5, +-1e300, +-2^-1074, +-1075, 1024, +-inf, +-NaN)
checks every case C Annex F defines bitwise against libm.

**On Linux the libm side is glibc, not macOS libm.** The results below were
recorded against macOS libm. Only C Annex F special cases are checked
bitwise against libm; other listed special points (tiny |x|) are checked
against MPFR, and a libm that isn't correctly rounded there is printed as a
`note:` line, not a failure (glibc's tanh at 0x1.fffffffffffffp-28 is 1 ulp
off; the recipe is correctly rounded, Codex E2E run 2026-10-19). The gate is the MPFR
column (at most 1 ulp, as below), the Annex F special cases, and the
`-O0`/`-O2` hash; compare any libm-only mismatch against MPFR before
blaming the recipe, and against the GPU bit check (`ocml_gpu/`), which is
the authority for bit-exactness.

## Results (Apple clang 21, arm64 macOS, 2,000,000 samples per sweep)

| Sweep | max ulp vs CR | 0 ulp | 1 ulp | libm max ulp | special mismatches |
|---|---|---|---|---|---|
| pow special grid (696) | 1 | 97.1% | 2.9% | 1 | 0 of 546 C-defined |
| pow wide (x 2^-1074..2^1024, \|y\| 1e-3..1e3) | 1 | 85.9% | 14.1% | 1 | 0 |
| pow random bit patterns | 1 | 99.8% | 0.15% | 1 | 0 |
| pow near-1 base, \|y\| 1..1e16 | 1 | 94.9% | 5.1% | 1 | 0 |
| pow integer y in [-60, 60], signed x | 1 | 75.0% | 25.0% | 1 | 0 |
| pow half-integer y | 1 | 78.2% | 21.8% | 1 | 0 |
| pow T9 1e-3..1e6, y in {2/3, 1/3, -3/2, -2/3, 5/6, 0.53, 3.263, 9.227, 0.25, 1.9885, 0.1, 5/3, 5/12} | 1 | 76.4% | 23.6% | 1 | 0 |
| pow triple-alpha args (0.025/T9)^3.263, (T9/0.025)^9.227 | 1 | 74.4% | 25.6% | 1 | 0 |
| pow VODE x^(1/L), L = 1..13 | 1 | 78.4% | 21.6% | 1 | 0 |
| pow VODE x^L, L = 1..13 | 1 | 76.1% | 23.9% | 1 | 0 |
| pow 10^t, t in [-12, 15] | 1 | 74.6% | 25.4% | 1 | 0 |
| pow overflow/underflow edges (results near 2^+-1024, subnormal) | 1 | 88.8% | 11.2% | 1 | 0 |
| log10 special values (26) | - | - | - | - | 0 |
| log10 random positive bit patterns | 1 | 99.985% | 0.015% | 1 | 0 |
| log10 near 1 | 1 | 99.46% | 0.54% | 1 | 0 |
| log10 1e-5..1e16 | 1 | 99.72% | 0.28% | 1 | 0 |
| tanh special values (26) | - | - | - | - | 0 |
| tanh [-25, 25] | 1 | 99.80% | 0.20% | 1 | 0 |
| tanh small 2^-40..1 | 1 | 97.9% | 2.1% | 1 | 0 |
| tanh random bit patterns | 1 | 99.96% | 0.04% | 1 | 0 |

macOS libm is correctly rounded in 99.6-99.999% of these samples. The
recipes never exceed 1 ulp; pow is faithful but not correctly rounded in
about a quarter of typical cases, which is OCML's design (`expep` rounds
`exp(x.hi)` before applying the `x.lo` correction). `log10(10^k)` is exactly
`k` for all `|k| <= 22`. `-O0` and `-O2` builds give identical bits.

Notable OCML behaviour reproduced by the recipe:
`pow(x, 1.0)` is not always `x` (6.2% of random `x`; for example
`pow(3, 1) = 0x1.7ffffffffffffp+1`), because OCML has no `y == 1` special
case. VODE evaluates `pow(x, 1.0 / L)` and `pow(x, L)` with `L = 1`.

## GPU bit check (`ocml_gpu/`)

Done on moth (gfx90a, 2026-10-18): the Loom recipes for pow, log10, tanh,
exp, log and cbrt are bit-identical to hipcc's OCML in every sweep (~1M
samples each), including `pow(x, 1) != x` for the same 6.21% of samples on
both sides and the `MATH_FAST_RCP` substitution. The only differences were 3
cbrt NaN payloads, which pass under the NaN rule. Stock `std::` calls and
direct `__ocml_*` calls gave identical bits everywhere. Rerun it when a
recipe, the OCML revision, or the GPU target changes.
