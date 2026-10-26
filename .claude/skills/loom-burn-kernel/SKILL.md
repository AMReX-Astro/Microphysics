---
name: loom-burn-kernel
description: Port a Microphysics burn (a chosen network + EOS + ODE integrator + screening/NSE/compile-time options) to a Loom (hrx-system) GPU kernel that is compiled with loom-import-cxx/loom-compile and launched from unit_test/test_react, then validate it against the stock HIP/CPU burner. Use this whenever the user wants to run a Microphysics network, integrator, or burner through Loom or HRX, make a "Loom kernel" or "Loom version" of a reaction network, compare Loom vs hipcc on test_react, or carry the chemistry-kernel autoresearch lessons over to Microphysics, even if they only name the network (e.g. "do aprox13 with VODE on Loom").
---

# Microphysics burn kernel → Loom kernel

The goal is a Loom kernel that does exactly what `do_react()` in
`unit_test/test_react/react_zones.H` does for one zone, for one fixed build
configuration, and a small opt-in change to `test_react` that loads the Loom
code object and launches it in place of the `amrex::ParReduce` burn.

This skill distills the chemistry-kernel work in hrx-system
(`experimental/loom-codegen-comparison/`, 12 autoresearch campaigns
`adv-mi210*`, 2026-10-10 → 10-18) and the first Microphysics port
(aprox13 + helmholtz + VODE, 2026-10-18). Almost every lost hour came from
process mistakes, not from the physics, so the gates below are deliberately
strict. Read `references/footguns.md` before touching anything; it is short
and each item cost real time once.

**Paths.** `<skill>` below means the directory holding this `SKILL.md`. In
Microphysics that is `.claude/skills/loom-burn-kernel` (Claude Code), also
reachable as `.agents/skills/loom-burn-kernel` (a symlink, for Codex); any
other install location works the same. Paths without `<skill>/` (e.g.
`references/footguns.md`) are relative to it; repo paths start at the
Microphysics root (`unit_test/...`).

**Agent environment (Claude Code or Codex).** Setup needs network access
(apt, git clone, pip), `sudo`, and detached jobs that run 15–30 min
(`nohup setsid ... &`, checked later). If the agent's sandbox blocks any of
those, have the user run the Setup block themselves, or switch to a mode
that allows them, before Phase 0. The `assets/prompts/` files are task
specs: with sub-agents, give one to each; without them, work through them
in order yourself.

**Never estimate the time.** Run `date -u` before every timestamp and every
time-based decision; agents in both end-to-end runs wrote guessed times into
their logs.

## What "done" means

1. Kernel source that compiles two ways from one file: under Loom
   (`__loom__` defined) and with the host C++ compiler (native path). The
   native path is how you prove the translation is right without a GPU.
2. Native differential validation: the native kernel matches the stock
   Microphysics burner zone by zone on the test_react grid, bit for bit.
3. Import + compile on the GPU host (this machine) with
   frozen, md5-recorded tools.
4. GPU correctness inside test_react: Loom vs the stock HIP build of the same
   configuration, same GPU, same inputs, bit for bit.
5. A hand-off note in `agent-notes/` (see Phase 8).

Performance work comes after all of that and is optional (Phase 7). A fast kernel that is not proven correct was
the single most expensive failure mode of the campaigns.

## Fixed policies (user decisions, 2026-10-18; don't ask about them)

- **Neutrinos off first, always.** The first milestone of every port is built
  with `USE_NEUTRINOS=FALSE` on the make command line (most networks'
  `Make.package` set it `TRUE`; the command line overrides it, which drops
  `-DNEUTRINOS` and the `neutrinos/` sources). Neutrino losses are a later
  milestone on top of a passing neutrinos-off kernel.
- **No FP contraction in reference builds.** Every build used as a correctness
  reference (stock CPU, stock HIP, the native path of the kernel, and any
  hipcc build of the kernel source) adds `XTRA_CXXFLAGS=-ffp-contract=off`.
  hipcc and GCC fuse multiply-adds by default, so without it the references
  don't compute what the source says and won't agree bit for bit. FMA may be
  enabled only in Phase 7, as an explicit optional experiment.
- **Bit-exact is the target**, on the host and on the GPU. Fall back to
  "close to machine precision" only where reassociation makes bit-exact
  impossible, and only with the user's approval and per-zone evidence. Any
  non-statement-for-statement rewrite needs that approval first.
- **NaN sign and payload bits don't matter**: any NaN matches any NaN in
  every comparison (`references/test-react-integration.md`, "Comparing").
- **`network.use_tables=1` is rejected at launch**; the prad limiter
  (`eos.prad_limiter_*`) is ported and supported.
- **No network-specific files in the skill's `assets/`**, not even as
  examples; ports live under `unit_test/test_react/loom/`.

## Setup — a fresh machine (once)

The skill assumes nothing is installed beyond ROCm. On a fresh Ubuntu 24.04
VM with ROCm. Everything runs on the GPU host itself; remote-host tooling is
in the `loom-burn-kernel-tools` follow-up. Tested end to end on 2026-10-18 for
`triple_alpha_plus_cago` (agent-notes `03-loom-e2e-container`):

```bash
sudo apt install build-essential git python3 python3.12-venv patch cmake ninja-build \
     g++-14 libstdc++-14-dev clang-20 clang-tools-20 llvm-20-dev libclang-20-dev
sudo usermod -aG render,video "$(id -un)"     # GPU access (log in again); agent shells may lack $USER
export ROCM_PATH=${ROCM_PATH:-/opt/rocm}
export LLVM_BIN=/usr/lib/llvm-20/bin            # inventory.py, CXX=clang++
export PATH=$ROCM_PATH/bin:$ROCM_PATH/llvm/bin:$LLVM_BIN:$PATH
GFX=$(amdgpu-arch | head -1); echo "${GFX:-no GPU: set GFX=gfx90a (or your target) by hand}"
export SKILL=$(realpath <skill>)                # absolute: the block below cd's away from the repo
export LOOM_WORK=${LOOM_WORK:-$HOME/loom-work}  # every script's base dir
mkdir -p "$LOOM_WORK" && cd "$LOOM_WORK"
python3 -m venv venv && venv/bin/pip install ast-grep-cli && export PATH=$PWD/venv/bin:$PATH
git clone https://github.com/AMReX-Codes/amrex.git && export AMREX_HOME=$PWD/amrex
git clone https://github.com/ROCm/hrx-system.git
git -C hrx-system checkout fa8eba4554           # pinned compiler (Phase 1); not 1468c06034+ (#2337)
(cd hrx-system && python3 dev.py cmake setup)   # .venv with cmake/ninja; see below if it exits 1
SRC=$PWD/hrx-system OUT=$PWD/loom-tools/pinned TARGETS=$GFX \
  nohup setsid bash $SKILL/scripts/build_loom_tools.sh > tools-build.log 2>&1 < /dev/null &
```

Notes from the end-to-end run:
- `uv` is optional: every Python script is standard library only, so
  `python3 <script>` works wherever this file says `uv run --no-project`.
  Ubuntu has no uv, pipx or npm by default, hence the venv for ast-grep.
- Ubuntu's `clang-20` installs `clang++` only under `$LLVM_BIN`, so keep it on
  `PATH` (`gen_network_data.py` runs `CXX=clang++`).
- `dev.py cmake setup` may exit 1 in a semgrep step (`semgrep-core: unknown
  option '--disable-version-check'`); harmless if `hrx-system/.venv/bin/cmake`
  and `.venv/bin/ninja` exist (all three end-to-end runs). It is the
  supported way to get hrx-system's pinned build tools; `build_loom_tools.sh`
  also accepts the apt `cmake`/`ninja-build` on `PATH` when there is no
  `.venv`, but that path has not been tested end to end.
- No GPU on the box (a build-only VM): Phase 6 stops after the production
  root compile, the HIP build and the FMA audit (`loom_gpu_compare.sh` skips
  its compare runs when no GPU is visible, or with `COMPILE_ONLY=1`).
  `amdgpu-arch` prints nothing. Set
  `GFX` by hand and pass it everywhere a target is asked for (`TARGETS`,
  `LOOM_TARGET`, `TARGET`, `ARCH`, `AMD_ARCH`); the scripts warn when they
  fall back to gfx90a.
- `g++-14` becomes the GCC that every clang (including hipcc) uses. Stock
  HIP builds of networks without intermediate rates then need the libstdc++
  workaround in `references/footguns.md` ("std::array<T, 0>") until
  Microphysics PR #2116 is merged.

Microphysics needs `AMREX_HOME` for every `make`, including the inventory's
`make -n`. Wire test_react once per checkout:
`git apply <skill>/assets/harness/test_react.patch`
(`references/test-react-integration.md`). The tool build checks for a C++23
libstdc++ (`std::expected`) first and takes 15–30 min (13 min at `-j8` under
Rosetta); start it first.

**Work order.** The phases are numbered by dependency, but the Loom tools
are the long pole, so: start the tool build (Setup), then do Phase 0 → 2 →
3 → 4 (all native, no Loom tools) while it runs. As soon as the tools exist,
run Phase 1 (probes), then smoke-compile each module once its own native
validator passes (EOS, then rates, then network; Phase 5 and NEW_NETWORK
step 7), even while later validators still run; the full-burn smoke and the
production root wait until the whole native suite has passed. If a probe
contradicts a translation choice, fix the translation before the smoke
compiles.

## Phase 0 — pin the configuration (no code yet)

Ask the user for anything else not given; do not guess these:
- `NETWORK_DIR`, `EOS_DIR`, `INTEGRATOR_DIR`, `SCREEN_METHOD`,
  NSE (`USE_NSE_TABLE`/`USE_NSE_NET`), and which `inputs_*` file.
  test_react is Strang (not SDC); say so if the user asks for SDC.
- GPU target: detect it (`amdgpu-arch`) and use it for both `AMD_ARCH` and
  `loom-compile --target`. `gfx90a` (MI210) is the proven target and
  `gfx942` (MI300X) was also used; anything else (e.g. wave32 RDNA) is
  unproven, so tell the user, rerun the OCML GPU bit check (in
  the `loom-burn-kernel-tools` follow-up) and watch per-lane scratch (~213 KB for
  aprox13). CDNA1/gfx908 is out of scope.
- The GPU must be on this machine; if it isn't, tell the user (remote-host
  tooling is in the `loom-burn-kernel-tools` follow-up).

Then record, in a new `$NOTES_DIR.md` (default
`NOTES_DIR=agent-notes/<NN>-loom-<network>-<integrator>`; set it to keep
notes elsewhere, and use `$NOTES_DIR/` wherever this file says
`agent-notes/<NN>-loom-<network>-<integrator>/`):
- Microphysics commit, AMReX commit, hrx-system commit
  (`git -C <hrx> rev-parse HEAD`), ROCm version (`cat /opt/rocm/core-*/.info/version
  /opt/rocm/.info/version 2>/dev/null; hipconfig --version`), GPU target, the
  exact make command line (including `USE_NEUTRINOS=FALSE` and
  `XTRA_CXXFLAGS=-ffp-contract=off`).
- The fully resolved compile flags. Get them from the real build, not from
  reading GNUmakefiles: `make -n <vars> | grep 'main.cpp'` (or build once
  with the vars and read the command line). Keep the `-D` list.
- The runtime parameters that the burn path reads and their values from the
  inputs file plus `_parameters` defaults. These become kernel inputs.

Work on a branch off `development` (AGENTS.md); never commit on
`development`, and do not commit at all unless the user asks.

## Phase 1 — probe the Loom frontend at the pinned commit

Loom's C++ importer changes weekly (e.g. #1993 added non-escaping value
references after the campaigns had rewritten every reference away). Do not
trust any list of "unsupported" constructs, including the one in
`references/loom-frontend-rules.md`, until you have re-probed it at the
commit you are using. Copy the skill's `probes/` kit (124 one-construct
probes with root `probe::k`, and `probe.sh`, which runs them all with one
frozen tool set) into the port's agent-notes directory, add probes for
anything the Phase 2 inventory found that the kit doesn't cover, and run it
on the GPU host as its header describes. `probes/README.md` has the measured
results at public hrx-system main `fa8eba4554` (2026-10-20, gfx90a), the
commit to check out (no patch is required). Member functions, `operator()`
accessors, `const` struct references, array reference parameters and deep
ternaries import, so translations don't rewrite them, and the shipped modules
use them (`Array1D`/`Array2D` accessors with the stock indices, `const&` for
everything read-only, pointers only for structs a function writes).
Uninitialized scalar and array locals, `for (;;)`, namespace-scope
`constexpr` tables (`aion[]`, rate tables) and nested multi-exit loops work
too. If you
use `fa8eba4554` on gfx90a, cite the table; on another commit or target,
rerun the kit and update it. Record the
tool md5 and hrx-system commit in the agent note.

The kit includes the f64 math probes (`m_f64_*`; Phase 2 tells you which
functions the burn needs) and a contraction probe (`m_fp_contract`: with
`--fp-model=precise`, `a * b + c` must compile to a separate multiply and
add, not `v_fma_f64`). Whatever math does not lower natively for f64 needs a
source recipe (`references/microphysics-translation.md`, "Math").

Build the frozen tools with `scripts/build_loom_tools.sh` (writes
`tools.md5`); no local patch is required at `fa8eba4554`.

## Phase 2 — map the reachable code

Run the inventory script first (it needs LLVM's `clang++` and `clang-query`:
Ubuntu `clang-20 clang-tools-20` (Setup) or Homebrew `llvm`, found via `--llvm` /
`LLVM_BIN`, `PATH` or `/usr/lib/llvm-*/bin`; plus `AMREX_HOME`. It runs on
CPU, takes 1–5 min, and runs no Loom tools). Its math list also tells Phase 1 which functions to probe.

```bash
uv run --no-project <skill>/scripts/inventory.py \
  --out agent-notes/<NN>-loom-<network>-<integrator>/inventory \
  -- NETWORK_DIR=... EOS_DIR=... INTEGRATOR_DIR=... USE_NEUTRINOS=FALSE
```

It builds clang's call graph and the closure from `do_react`
(`reachable.txt`) with `main.cpp`'s flags from `make -n`, and runs
clang-query matchers for the constructs in `references/loom-frontend-rules.md`.
Output: `summary.md` (counts, globals that become `BurnParams`/tables,
`constexpr` arrays, math calls, "CPU vs HIP defines", per-function table) and
`hits.tsv` (every hit with file:line, function, reachable flag). Its limits
(reachability is an upper bound, runtime-gated paths, CPU vs HIP defines,
what it doesn't measure) are in `references/microphysics-translation.md`,
"Inventory limits"; cover them by reading code.

Put the summary and your reading of it in the agent note: the reachable
functions by module, the hardest constructs, the parameters and tables,
the math functions needing recipes, and a size estimate in translated lines.
If that is far larger than the chemistry kernel (~3k translated lines), tell
the user before starting. Autodiff screening and rates are the main size
driver left once neutrinos are off; even aprox13 + helmholtz + VODE without
neutrinos is several times the chemistry kernel.

## Phase 3 — translate

**For an `rhs.H`-framework network with helmholtz + VODE, follow
`assets/NEW_NETWORK.md` step by step**: it says what is reused as is
(`assets/kernel/`, `assets/modules/` including the generic `network_rhs.H`),
what is generated (`network_data.H`, `burn.H` via
`scripts/gen_network_data.py`, see `references/microphysics-translation.md`,
"Generated network data"), and what is hand-translated per network
(`rates.H`, `network.H`, with the `assets/prompts/` templates). For anything
else (another EOS, integrator, NSE, pynucastro networks) read
`references/microphysics-translation.md` and
`references/loom-frontend-rules.md` first.

Work from the source plus `dual.H`'s guide and the differential validators.
For a large network (dozens of rates, deep templates), C++ Insights'
specialized source helps (`specialize.py`, in the `loom-burn-kernel-tools` follow-up).

Lint every translated file before any Loom import (local, well
under a second; rules and caveats in `references/loom-frontend-rules.md`,
"Lint"):

```bash
K=unit_test/test_react/loom/<network>_<integrator>
ast-grep scan -c <skill>/scripts/loom-lint/sgconfig.yml \
  $(ls $K/*.H $K/*.cpp | grep -v '/validate_')
```

Never lint the `validate_*.H` files: they are host-only and use std,
lambdas and references on purpose (hundreds of hits). The shipped
`network_rhs.H` gives 9 `warning[no-reference]` hits (the `rate_t&` outputs
of `fill_additional_rates`, written only by whole-struct assignment, which
imports); those are expected. Your own files must have no errors, and every
`no-reference` warning in them must be a reference the probes showed
import.

Core rules:

- **One dual-compile source tree** in
  `unit_test/test_react/loom/<network>_<integrator>/`, started from the
  skill's `assets/` (see "Assets" below and `assets/README.md`).
  `#if defined(__loom__)` selects Loom bindings; otherwise host math.
- **No AMReX or std headers in the Loom path.** `using Real = double;` inside
  the kernel namespace stands in for `amrex::Real`; the host-side code in
  test_react keeps using `amrex::Real`.
- **Preserve the computation**: same species, rates, screening, EOS,
  tolerances, step-size control, iteration limits, and floating-point
  evaluation order. Translate statement for statement. Don't specialize on
  the test grid (no hard-coded densities, temperatures, cell counts).
- **Prefer a generator over hand edits** when the source is generated
  (pynucastro networks) or large: a script that extracts bodies from the
  pinned Microphysics source, applies exact-match replacements, and fails if
  the source digest changes.
- **Kernel ABI**: pointers plus scalars only; no by-value records and no
  pointer to a struct of pointers. Runtime parameters arrive in one
  `const BurnParams*`; each table is its own pointer argument. Mirror every
  record that crosses the boundary on both sides and `static_assert` its
  `sizeof`/`offsetof` in both the host and kernel builds:
  `scripts/gen_abi_layout.py` measures every `*Params` struct on the host and
  writes `abi_layout.H`, which the shipped kernel root includes (the lint
  allows `__builtin_offsetof`).

## Phase 4 — native validation (local, cheap, do it before any GPU work)

Build the native path into a CPU test_react (the `assets/harness/` hook and
`assets/validators/`; `USE_LOOM_VALIDATE=TRUE LOOM_KERNEL_DIR=loom/<port>`) and compare against the stock burner zone by zone,
built with `XTRA_CXXFLAGS=-ffp-contract=off` on both sides. The CPU
reference must take the same physics path as the HIP build: neutralize every
physics define in the inventory's "CPU vs HIP defines" (for VODE,
`ALLOW_JACOBIAN_CACHING`: see `references/test-react-integration.md`,
"Native validation"). Compare over the full test_react grid (`n_cell` ≥ 16 so
the grid spans the whole T/ρ/X range), double bit patterns and integer
counters (`n_rhs`, `n_step`, success). Also validate subroutines in
isolation (RHS, Jacobian, EOS call, LU/solve, rates, duals); a whole-burn
mismatch is hard to localize. Validate in translation order: the EOS and
dual validators need only the copied modules, the generated
`network_data.H` and `harness/loom_eos_params.H`, so they pass before any
rates or network code is written (`assets/NEW_NETWORK.md`, step 3). Zero mismatches is the bar. The f64 math
recipes use host libm on the native path, so they are validated separately
(host accuracy and GPU bit check, in the `loom-burn-kernel-tools` follow-up).

## Phase 5 — import and compile on the GPU host

Only send sources that pass the Phase 3 lint. Run full-kernel
`loom-import-cxx`/`loom-compile` only on the GPU host, one at a time: a full
burn needs 8+ GB, and parallel compiles on a Mac caused a kernel panic.
Drive every import/compile through the scripts:
`scripts/loom_smoke.sh launch|status|wait|fetch` smoke-compiles each
module as soon as it validates (templates in `assets/smoke/`; `wait` keeps
compiles one at a time), and `scripts/loom_gpu_compare.sh` compiles the
production root. Commands and flags are in
`references/test-react-integration.md`. Keep the compile report (the scripts
write it to `<root>.report.json`, not the log) and VGPR/AGPR/scratch
from every build (`llvm-readelf --notes`; the scripts print them). A fused full burn
compiles in 9 min (triple_alpha) to about 30 min (aprox13), almost all in
register repair:
use tools built from the pinned commit and start the production root compile as
soon as the burn body compiles (`references/footguns.md`).

## Phase 6 — GPU correctness in test_react

Apply the host-side change described in `references/test-react-integration.md`
(opt-in via a runtime parameter; the stock path stays the default; start from
`assets/harness/loom_launch.H` and the `loom_params.H` template). Build
test_react with `USE_HIP=TRUE AMD_ARCH=<target> USE_LOOM_BURN=TRUE
LOOM_KERNEL_DIR=loom/<port> USE_NEUTRINOS=FALSE XTRA_CXXFLAGS=-ffp-contract=off`
(`scripts/build_test_react.sh loom/<port>` runs it detached;
`build_test_react.sh status` shows it, and the log
is `$LOOM_WORK/test_react-build/last/build.log`, with the hipcc
version and any `HIPCC_*_FLAGS_APPEND` at the top). For a network without intermediate rates on
libstdc++ >= 13 (Ubuntu 24.04), the stock HIP build needs the
`std::array<T, 0>` workaround in `references/footguns.md` until Microphysics
PR #2116 is merged.

Check for FP contraction with `scripts/fma_audit.py --loom loom_root.loom
--hsaco loom_root.hsaco`. A raw `v_fma_f64` count is not evidence: AMDGPU
lowers every f64 division to 5 `v_fma_f64` (with `v_div_scale`/`v_div_fmas`/
`v_div_fixup`) and every f64 sqrt to a few more, on top of the recipes'
explicit `fma` calls (triple_alpha: 32,866 total = 19,401 explicit + 5 x
2,495 divisions + ~990 for 303 sqrts). The script does that accounting as a
screening heuristic, not an exact attribution: it allows up to 8 remaining
FMAs per sqrt and treats a negative remainder as constant folding, so a few
contracted operations could hide in that slack. Report its result as "no
sign of contraction", and rest the claim on the `m_fp_contract` probe
(Phase 1: separate `v_mul_f64` + `v_add_f64`) together with it. If it reports
unexplained FMAs, stop and tell the user. Then, on the same GPU in the same process
(`scripts/loom_gpu_compare.sh`):
1. run the stock burn, then the Loom burn, at least twice each, and check
   that each side repeats bit-identically;
2. compare every output component zone by zone, bit for bit (NaN rule
   aside); report the worst zone and per-component differences on failure;
3. also compare `n_rhs`/`n_step`; a different step count with matching
   state is a sign of a changed control path, not noise;
4. include a run with every supported non-default path the inputs don't
   exercise (e.g. the prad limiter on: `loom_gpu_compare.sh` with
   `RUNS=";eos.prad_limiter_rho_c=1.0e7 eos.prad_limiter_delta_rho=3.0e6"`).

If something fails, localize it first (native path? recipe math? ABI?
layout?). Accepting a non-bit-exact result is the user's decision, with the
evidence in front of them.

**Standalone HIP reproducer.** The kernel tree also compiles with hipcc as an
ordinary HIP kernel (`-DMPL_HIP_STANDALONE`, `support.H`). The case-file
writer for test_react, the driver that burns a case file outside
Microphysics and the vendor-reproducer bundler are in the `loom-burn-kernel-tools` follow-up.

## Phase 7 — performance (optional, only after Phase 6 passes)

Timing must be paired and interleaved against the reference on the same
GPU, with clock/temperature telemetry; a single "Run time" line from
test_react is not evidence. FMA contraction is tried only here. The
performance guide and the compiler-diagnostic tools are in
the `loom-burn-kernel-tools` follow-up.

## Phase 8 — hand-off

Update the agent note with: configuration and commits, tool md5s, probe
results, every rewrite that is not statement-for-statement and why, every
non-bit-exact change and who approved it, validation commands and results,
GPU results, compile stats, open bugs with Loom issue numbers, and the next
steps. If you found a Loom bug, write a minimal reproducer and offer to file
it; don't file it without asking. When a port improves a reusable asset,
copy it back into `assets/` (network-neutral only) and update its README.
Check the time with `date` before writing any timestamp.

## References

- `references/footguns.md` — read first. Process, machine and compiler hazards (incl. compile time / register repair).
- `references/loom-frontend-rules.md` — construct → rewrite table as of the last probes, `MPL_ASSUME`/SUBRANGE, lint, kernel root shape.
- `references/microphysics-translation.md` — inventory limits, burn_t, runtime params and unported paths, tables, generated network data, autodiff, math, NSE.
- `references/test-react-integration.md` — layout, kernel signature, host launcher, native validation, build/run and import/compile commands, comparison and NaN rule.
- `probes/README.md` — Phase 1 probe kit (`probe.sh`, `src/`) and results at `fa8eba4554`.

## Scripts

- `scripts/inventory.py` (+ `reachable.py`) — Phase 2 call-graph closure and construct inventory.
- `scripts/gen_network_data.py` + `extract_network_data.cpp` — generate `network_data.H` / `burn.H` for `rhs.H`-framework networks, checked against the stock source.
- `scripts/loom-lint/` — ast-grep rules for Loom-hostile constructs (`tests/selftest.sh`).
- `scripts/build_loom_tools.sh` — build frozen `loom-import-cxx`/`loom-compile` from an hrx-system tree plus optional patches; writes `tools.md5`. Header has the fresh-machine prerequisites.
- `scripts/loom_smoke.sh` — launch/status/wait/fetch a Loom import + compile of a smoke kernel on the GPU host. `wait` serializes compiles; `wait`/`status` exit 6 when an import or compile failed.
- `scripts/gen_abi_layout.py` — Phase 3 ABI checks: host-measured `sizeof`/`offsetof` of every `*Params` record, written to the port's `abi_layout.H` (included by the kernel root).
- `scripts/fma_audit.py` — Phase 6 contraction screen: explains the `v_fma_f64` count (explicit recipe FMAs from the `.loom` IR, 5 per f64 division, a bounded remainder per sqrt); a heuristic, paired with the `m_fp_contract` probe.
- `scripts/build_test_react.sh` — detached HIP build of test_react with the Loom launcher, in place; `status` shows the log tail.
- `scripts/loom_gpu_compare.sh` — Phase 6 on the GPU host (run it directly when local): production root compile, kernel-argument metadata, test_react compare runs.

## Assets (`assets/README.md`)

- `assets/NEW_NETWORK.md` — step-by-step bring-up of another `rhs.H`-framework network.
- `assets/kernel/` — platform layer `support.H` (`MPL_DEVICE`, `mpl::math`, `MPL_ASSUME`), OCML f64 recipes `f64_math.H` / `f64_math_extra.H`, POD autodiff `dual.H`.
- `assets/modules/` — network-independent helmholtz `eos.H`, `screen.H`, `network_rhs.H`, `vode.H`, `burner.H`.
- `assets/root/burn_kernel.cpp` — kernel root `mpl::burn_zones`.
- `assets/harness/` — `test_react.patch` (test_react wiring), native validation hook, dual validator, `loom_launch.H` launcher, `loom_params.H` template (host `fill_*` and unsupported-setting checks).
- `assets/validators/` — network-neutral validator templates (EOS, rates, network, whole-burn integrator, kernel root) and `validate_port.H` (per-port hooks).
- `assets/smoke/` — smoke-compile kernels for `loom_smoke.sh`.
- `assets/prompts/` — agent prompt templates (`_common.md` + per task).

## Moved to the follow-up

The `loom-burn-kernel-tools` follow-up holds the optional tooling:
remote-GPU-host support (`LOOM_HOST`, `ship_test_react.sh`), the standalone
HIP reproducer and its test_react case-file writer (`make_standalone.sh`,
`loom_case.H`), the f64 recipe accuracy tests, the performance guide,
compiler diagnostics (`trace-repair.patch`, `repair_summary.py`,
`loom_stack_sampler.sh`, `compare_loom_tools.sh`, `hsaco_stats.py`),
`specialize.py` (C++ Insights) and the seven-dwarfs test configurations.
