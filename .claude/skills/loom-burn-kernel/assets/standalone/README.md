# Standalone burn reproducer (HIP, no Microphysics or AMReX)

This directory is a self-contained GPU program that runs one Microphysics
reaction-network burn: a nuclear network, the helmholtz equation of state and
the VODE stiff ODE integrator, translated statement for statement from
[Microphysics](https://github.com/AMReX-Astro/Microphysics) into plain C++
(the "kernel tree": `burn_kernel.cpp` and the `.H` files next to it). It
needs only ROCm's `hipcc`. It was made to compare compilers on the same
source and to hand a complete, small reproducer to a compiler or GPU vendor.

The kernel tree compiles three ways from one source:

| Build | Selected by | Used for |
|---|---|---|
| hipcc, this directory | `-DMPL_HIP_STANDALONE` | this program |
| Loom (hrx-system) | `__loom__` (set by `loom-import-cxx`) | the Loom code object (`--hsaco`) |
| host C++ compiler | neither | bit-for-bit validation against Microphysics on the CPU |

A **case file** (`cases/*.case`) holds everything the burn reads (runtime
parameters, the 14 helmholtz EOS tables, the initial state of a box of
zones) and a reference result from a known-good build (usually the stock
Microphysics HIP build on the same GPU). Its layout is documented in
`loom_case_format.H`.

Use a case written by a **GPU** build to check a GPU binary. test_react
builds its initial ρ/T grid with `pow` and `log10` in a device loop, so a
CPU-written case has slightly different *inputs* (glibc vs OCML rounding)
as well as a CPU reference: the standalone kernel still burns those inputs
faithfully, but it is not expected to match the CPU reference bit for bit.

## Build

```bash
make                    # burn_standalone_ocml: hipcc's own device math library (OCML)
make MATH=recipes       # burn_standalone_recipes: the f64 math transcriptions the Loom build uses
make ARCH=gfx942        # target (default: the first GPU amdgpu-arch reports, else gfx90a)
make asm                # device assembly: burn_standalone_$(MATH).s
```

The flags are `-std=c++20 -O3 -ffp-contract=off`. Keep `-ffp-contract=off`:
the reference build uses it, and FMA contraction changes the bits.
`MATH=ocml` calls the same device math as the stock Microphysics HIP build;
`MATH=recipes` calls statement-for-statement transcriptions of OCML's f64
`exp`, `log`, `cbrt`, `pow`, `log10` and `tanh` (`f64_math.H`,
`f64_math_extra.H`) with `shim/loomcxx/scalar.h` mapping their primitive
operations to AMDGPU builtins. Both are expected to reproduce the reference
bit for bit.

## Run

```bash
./burn_standalone_ocml cases/<case>                     # hipcc kernel vs the case's reference
./burn_standalone_ocml cases/<case> --hsaco burn.hsaco  # also a Loom code object of the same root
./burn_standalone_ocml cases/<case> --zone 3,7,12       # one zone only (a minimal reproducer)
```

| Option | Meaning |
|---|---|
| `--repeats N` | burns per kernel (default 2); repeated burns must be bit-identical |
| `--no-hip` | skip the hipcc kernel (only with `--hsaco`) |
| `--hsaco FILE`, `--kernel NAME` | Loom code object and its symbol (default `mpl.burn_zones`) |
| `--loom-blocks N` | workgroups the code object was compiled for (default: zones of the case / 128, rounded up) |
| `--zone I,J,K` | burn only that zone of the case (global AMReX index, as printed in mismatch reports) |
| `--write-subset FILE` | with `--zone`: write that one-zone case (inputs + its reference) and exit |
| `--out FILE`, `--out-loom FILE` | write the case with the hipcc / Loom result as its reference |
| `--print-inputs` | print each burned zone's ρ, T and mass fractions as hex floats |

Every double is compared by its bit pattern (any NaN matches any NaN), and
`n_rhs`, `n_step` and the per-zone success flag must match exactly. Output:

```
RUN kind=hipcc run=0 ms=2512.3 failures=0
REPEAT kind=hipcc runs=2 identical=1
COMPARE hipcc_vs_reference a=reference b=hipcc zones=4096 identical_zones=4096 ... bit_exact=1
STANDALONE_RESULT PASS
```

On a mismatch each differing state component gets one line with the count of
differing zones, the largest relative difference, and the worst zone with
both values in decimal and hex, followed by the first differing zone and the
`--zone` option that reruns it alone. Exit status: 0 = everything compared
was bit-exact and repeatable, 1 = mismatch, 2 = error.

## Workflows

- **Is hipcc right?** `./burn_standalone_ocml cases/<case>` compares the
  hipcc-compiled kernel tree with the stock Microphysics HIP result. A
  mismatch is either a translation bug in the kernel tree (rule this out
  first: the tree matches Microphysics bit for bit on the CPU) or a
  compiler difference.
- **Loom vs hipcc on one GPU, one process:** `--hsaco` runs both on the
  same device buffers and prints `loom_vs_reference` and `loom_vs_hipcc`.
- **Minimal reproducer:** take the first differing zone from the report,
  `--zone I,J,K --write-subset one.case`, and ship this directory with
  `cases/one.case`. Someone else runs `make && ./burn_standalone_ocml
  cases/one.case`.
- **Compare two compilers without Microphysics:** run with compiler A and
  `--out a.case`; then build with compiler B (`make HIPCC=...`) and run it on
  `a.case`.

## Where this comes from

The bundle is made by the `loom-burn-kernel` agent skill in Microphysics
(`.claude/skills/loom-burn-kernel/scripts/make_standalone.sh`); `MANIFEST`
lists the source port, commits and file checksums. Case files are written by
Microphysics' `unit_test/test_react` built with the skill's harness, run with
`unit_test.loom_dump_case=<file>`.
