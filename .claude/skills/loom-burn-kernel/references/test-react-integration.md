# Running the Loom kernel inside unit_test/test_react

The change to test_react is small, opt-in, and leaves the default build and
its CI benchmarks untouched. The skill ships it as `assets/harness/`
(`test_react.patch` for GNUmakefile/Make.package/_parameters/main.cpp,
`loom_launch.H`, `loom_params.H` template, validation hook) and
`assets/root/burn_kernel.cpp`. In a fresh checkout, from the repo root:

```bash
git apply <skill>/assets/harness/test_react.patch
cp <skill>/assets/harness/{loom_launch.H,loom_validate.H,loom_validate.cpp,loom_eos_params.H} unit_test/test_react/
cp <skill>/assets/harness/loom_params.H unit_test/test_react/   # then fill its PORT: lines
```

 the port's own `loom/<cfg>/README.md` holds
the exact signature and commands for that port.

## Layout

```
unit_test/test_react/
  loom/<network>_<integrator>/   # one dual-compile kernel directory per port
    support.H, f64_math.H, f64_math_extra.H, dual.H   # from assets/kernel/
    eos.H, screen.H, network_rhs.H, vode.H, burner.H  # from assets/modules/ (BurnParams in burner.H)
    burn_kernel.cpp                                   # from assets/root/: root mpl::burn_zones
    abi_layout.H                 # generated (scripts/gen_abi_layout.py): params sizeof/offsetof
    network_data.H, burn.H       # generated (scripts/gen_network_data.py)
    rates.H, network.H           # hand-translated per network (assets/NEW_NETWORK.md)
    validate_*.H                 # from assets/validators/ + harness/validate_dual.H;
                                 # validate_port.H holds the per-network PORT hooks
    README.md
    tests/import/                # smoke kernels from assets/smoke/
  loom_eos_params.H    # as is: fill_eos_params / fill_eos_tables (all the EOS validator needs)
  loom_params.H        # per port, from the assets/harness template: host fill_* functions
                       # and loom_require_supported_settings(); shared by validators and launcher
  loom_launch.H        # host: load .hsaco, pack args, launch, compare
  loom_validate.H/.cpp # native validation hook (dispatches to the kernel dir's validate_*.H)
  GNUmakefile, Make.package, main.cpp, _parameters   # from test_react.patch

```

Records shared by host and kernel (`BurnParams` and its parts) use only
`double`, `int`, `long` and pointers, with `static_assert`s on
`sizeof`/`offsetof` on both sides.

## Kernel signature (pointers and scalars, no records by value)

Never pass a pointer to a struct that holds pointers (probe c06): each table
is its own pointer argument, and the kernel rebuilds the tables struct by
value (c04). A pointer to a struct of scalars (`const BurnParams*`) is fine
(c01). The shipped root (`assets/root/burn_kernel.cpp`) takes `params`, the
14 helmholtz tables, the state box (`Array4` pointer, `jstride`, `kstride`,
`nstride`, box lower corner and extents), the `plot_t` component indices, the
`n_rhs` box (2 components: `n_rhs`, `n_step`) with its strides, `dt`, and a
per-zone `zone_success` array that the host reduces into `num_failed` (no
atomics). One lane does one zone and reproduces `do_react()` exactly: read
ρ, T, X, `normalize_abundances_burn`, `e = 0`, `T_fixed = -1`, `burner`, then
write `ispec`, `irodot`, `irho_hnuc` and `n_rhs` comps 0/1. The host form of
the root takes the zone index as a trailing `int`; the launcher
`static_assert`s its argument variables against it.

The ParReduce in `main.cpp` also returns the zone with the most RHS calls
(written to `zone_state.out`). The launcher finds that zone on the host from
`n_rhs` and re-burns it with the stock `do_react` on the host, so
`zone_state.out` holds the stock CPU result.

## Host side (`loom_launch.H`)

Compiled only when `USE_HIP=TRUE USE_LOOM_BURN=TRUE` (which adds
`-DMICROPHYSICS_LOOM`); `LOOM_KERNEL_DIR=loom/<port>` on the make line
becomes `-DLOOM_KERNEL_DIR_PATH=loom/<port>`, unquoted because the build
scripts receive `$(DEFINES)` inside `"..."`; the harness headers stringize
it (`LOOM_KERNEL_PATH(burner.H)`, in `loom_eos_params.H`). The patched
GNUmakefile stops with an error without it, or if the path contains
`linux`/`unix` (predefined macros).

0. `loom_require_supported_settings()` (in `loom_params.H`) aborts on
   runtime-gated paths the kernel does not implement
   (`microphysics-translation.md`, "Runtime-gated paths"), e.g.
   `network.use_tables=1`.
1. After `network_init()`, fill `BurnParams` from the `*_rp::` namespaces and
   copy it to an `amrex::The_Arena()` allocation. The helmholtz tables are
   `AMREX_GPU_MANAGED`, so their addresses are passed directly (worked on
   an MI210 with XNACK off). The launcher aborts if
   `hipDeviceAttributeManagedMemory` is 0; a device-copy fallback is not
   implemented (the stock HIP burn needs managed memory too).
2. `hipModuleLoad` the code object and `hipModuleGetFunction` with
   `unit_test.loom_kernel` (default `mpl.burn_zones`). Read the symbol from
   `llvm-readelf --notes` / `--syms` on the code object; do not guess it.
3. Build the argument array (`void* args[] = {&p0, &p1, ...}`, each scalar in
   a variable of the exact declared width) and `hipModuleLaunchKernel(fn,
   blocks,1,1, 128,1,1, 0, stream, args, nullptr)`. `blocks` must equal the
   compiled `workgroup_count.x`: the box's zone count is checked against
   `unit_test.loom_ncells`, so run with `max_grid_size = n_cell` (one box).
4. Use AMReX's current stream (`amrex::Gpu::gpuStream()`) and
   `amrex::Gpu::streamSynchronize()` so ordering with AMReX kernels holds.
5. Unload the module at the end.

## main.cpp change (in `test_react.patch`)

```cpp
    ValLocPair<int, burn_t> r;

#ifdef MICROPHYSICS_LOOM
    if (!unit_test_rp::loom_hsaco.empty()) {
        BL_PROFILE("do_react_loom");
        r = loom_burn_main(state, integrator_n_rhs, vars, names, num_failed_d);
    } else
#endif
    {
        BL_PROFILE("do_react");
        /* existing ParReduce, unchanged */
    }
```

plus `#include <loom_launch.H>` under `MICROPHYSICS_LOOM`, and, under
`LOOM_VALIDATE`, `#include <loom_validate.H>` and a hook right after
`network_init()` that runs `loom_validate()`, aborts on mismatches and
returns when `unit_test.loom_validate` is set.

New `_parameters` entries (namespace `unit_test`): `loom_validate` (int,
run the native validators instead of the burn), `loom_hsaco` (string,
empty = stock path), `loom_kernel` (symbol), `loom_ncells` (compiled zone
count), `loom_compare` (int, 0/1), `loom_repeats` (int).

## Native validation (Phase 4)

Build the same tree without `USE_HIP`, with the same `USE_NEUTRINOS=FALSE
XTRA_CXXFLAGS=-ffp-contract=off`, and `USE_LOOM_VALIDATE=TRUE
LOOM_KERNEL_DIR=loom/<port>` plus `XTRA_CPPFLAGS="-DLOOM_VALIDATE_<MODULE> ..."`
for the validators to run (`assets/README.md`, "harness/"):

```bash
cd unit_test/test_react
make -j NETWORK_DIR=... EOS_DIR=... INTEGRATOR_DIR=... USE_NEUTRINOS=FALSE \
     XTRA_CXXFLAGS=-ffp-contract=off USE_LOOM_VALIDATE=TRUE LOOM_KERNEL_DIR=loom/<port> \
     TMP_BUILD_DIR=tmp_val_<set> EBASE=loomval_<set> \
     XTRA_CPPFLAGS="-DLOOM_VALIDATE_DUAL -DLOOM_VALIDATE_EOS ..."
./loomval_<set>3d.gnu.ex inputs_<net> n_cell=16 integrator.use_jacobian_caching=0 unit_test.loom_validate=1
```

Use the default `COMP=gnu` for every native validation build (the same
compiler as the stock CPU reference and CI), and one `TMP_BUILD_DIR`/`EBASE`
per validator set: AMReX's make does not rebuild on flag changes. A fresh
`TMP_BUILD_DIR` can fail once with "No rule to make target
extern_parameters.cpp" under GNU make 3.81 (macOS); rerun make. The kernel
root is an ordinary inline function there, called zone by zone.

### Sensitivity check

Each validator prints its own sensitivity line (a 1-ulp perturbation that
must be detected). The manual check `prompts/_common.md` asks for is done in
a copy, so the real port is never touched (E2E run 2026-10-19):

```bash
cp -r loom/<port> loom/<port>_sens
# perturb ONE non-commutative operation in loom/<port>_sens, e.g. reassociate
# a product: a * (b * c)  ->  (a * b) * c   (a + b -> b + a is a no-op)
make -j NETWORK_DIR=... <same vars> LOOM_KERNEL_DIR=loom/<port>_sens \
     TMP_BUILD_DIR=tmp_val_sens EBASE=loomval_sens XTRA_CPPFLAGS="-DLOOM_VALIDATE_<MODULE>"
./loomval_sens3d.gnu.ex inputs_<net> unit_test.loom_validate=1 integrator.use_jacobian_caching=0
# expect LOOM_MISMATCH lines for the perturbed function; then
rm -rf loom/<port>_sens tmp_val_sens loomval_sens3d.gnu.ex
```

The CPU reference must take the HIP build's physics path. For VODE that means
`-DALLOW_JACOBIAN_CACHING` (CPU only): run both the CPU reference and the
native path with `integrator.use_jacobian_caching=0`, which makes `dvjac`
evaluate every time, the same as the HIP build where the caching code is
compiled out. (`XTRA_CPPFLAGS=-U...` does not work: it lands before the `-D`
on the compile line.) The kernel follows the HIP branch of every such
`#ifdef`.

## Build and run (GPU)

Stock reference and Loom launcher are the same binary:

```bash
cd unit_test/test_react
GFX=$(/opt/rocm/llvm/bin/amdgpu-arch | head -1)   # the same target for hipcc and loom-compile
make -j USE_HIP=TRUE AMD_ARCH=$GFX NETWORK_DIR=... EOS_DIR=... INTEGRATOR_DIR=... \
     USE_NEUTRINOS=FALSE XTRA_CXXFLAGS=-ffp-contract=off USE_LOOM_BURN=TRUE LOOM_KERNEL_DIR=loom/<port>
./main3d.hip.HIP.ex inputs_<net> n_cell=16 max_grid_size=16 integrator.use_jacobian_caching=0 \
     unit_test.loom_hsaco=burn.hsaco unit_test.loom_ncells=4096 unit_test.loom_compare=1
```

(Flag policy: SKILL.md, "Fixed policies".) `scripts/build_test_react.sh`
runs this build detached; `scripts/loom_gpu_compare.sh`
compiles the root and runs the compare runs.

Runtime library: by default run on ROCm's `libamdhip64` (AMReX's managed
memory works there). If you switch to HRX's (`ninja
libhrx_src_binding_hip_amdhip64`, pointed to via `LD_LIBRARY_PATH` and a
RUNPATH on the binary), first make AMReX run without managed memory
(`amrex.the_arena_is_managed=0`) and confirm the stock path still passes.
Record which runtime each result used.

## Import/compile commands

Run on the GPU host only, from a bash script, with the pinned Loom 0.1.0
tools (`hrx-loom==0.1.0`, or `scripts/build_loom_tools.sh` from the
`v0.1.0` tag; SKILL.md Setup):

```bash
loom-import-cxx --data-model=lp64 --fp-model=precise \
  --root=mpl::burn_zones --output=burn.loom loom/<cfg>/burn_kernel.cpp
loom-compile burn.loom --root=mpl.burn_zones \
  --config=mpl.burn_zones.workgroup_count.x=$(( (NCELLS + 127) / 128 )) \
  --config=mpl.burn_zones.workgroup_count.y=1 \
  --config=mpl.burn_zones.workgroup_count.z=1 \
  --format=amdgpu-hsaco --target=amdgpu:$GFX \
  --compile-report=summary --compile-report-output=burn.report.json --output=burn.hsaco
```

Wrap both in `timeout` and a CPU set (`taskset`) that avoids the cores
reserved for timing. The import takes seconds; a full burn's compile takes
about 30 min and 8+ GB, so give it a 3600 s timeout (why, and how to watch
it: `footguns.md`, "A full fused burn takes ~30 min"). (The kernel's
includes are all `"..."` next to the source, so no `--I` is needed.)
At the pin, aprox13 + helmholtz + VODE on gfx90a uses 256 VGPR and about 72 KB of
scratch per lane.

Keep from every build: the `.loom` IR (gzip), the compile report (without
`--compile-report-output` the JSON goes to stderr and makes the log ~84 MB;
gzip it), the disassembly, and VGPR/AGPR/scratch (`llvm-readelf --notes`).

## Comparing outputs

With `loom_compare=1` the launcher runs `max(2, loom_repeats)` stock/Loom
pairs from the same initial state into separate MultiFabs, checks that each
side repeats bit-identically, and compares stock vs Loom per component and
per zone: the bit pattern of every state component (`ispec`, `irodot`,
`irho_hnuc`, ...), `n_rhs`/`n_step`, and failure counts. It prints one
machine-readable `LOOM_COMPARE` line plus a per-component table (count of
differing zones, max abs/rel difference, worst zone with its ρ, T) and aborts
unless every zone is bit-identical, so scripts can gate on it. You can also
write both plotfiles and use AMReX `fcompare` as an independent check.

The target is bit-exact (SKILL.md, "Fixed policies"); there is no tolerance
knob. If bit-exact is impossible because of an approved reassociation, the
fallback ("close to machine precision", per-zone diffs reported) is the
user's decision.

**NaN rule** (user, 2026-10-18): NaN sign and payload bits don't matter; any
NaN matches any NaN. Implemented in `loom_same_bits` (`loom_launch.H`) and
`loom_check_bits` (`loom_validate.H`); the OCML GPU check treats NaN-payload
differences as passes. Every other value, including signed zeros and
infinities, must match bit for bit.
