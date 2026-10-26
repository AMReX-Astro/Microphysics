# Footguns from the chemistry-kernel campaigns

Each item happened at least once during the hrx-system chemistry work
(2026-10-10 → 10-18). The raw evidence (hrx-system `.notes/`, the campaign
branches) is in the original author's local trees only, not in any public
repository; the items below are the distilled result.

## Machine safety

- **No full-kernel Loom compiles on a laptop.** A Loom compile of a kernel this
  size takes several GB (8+ GB for a full burn); parallel local compiles
  panicked a Mac (2026-10-15). Forward-differentiating one large function made
  `loom-import-cxx` run 8+ minutes and was the likely memory hog. Import and
  compile the real kernel only on the GPU host: when that is the machine
  you are on (a Linux box or VM with ROCm and the GPU), use the scripts'
  `LOOM_HOST=local` mode (the default), check `free -g` first, and still run
  one full compile at a time. On a laptop without the GPU, tiny one-construct
  probes may run locally, one at a time, if the user agrees.
- **Remote host etiquette (moth; skip when the GPU is local).** Use only the plain `ssh moth` alias, one
  connection at a time, spaced out. Bursts of connections got sshd to refuse
  sessions for ~5 min and look like probing. Long work goes under `nohup`
  with a log; check on it with one long-lived session, not repeated
  reconnects. Never ask to bypass the firewall or jump through another host.
- **moth reboots for kernel updates** (e.g. 2026-10-18 08:47 AEDT, mid-compile;
  `last -x` shows the kernel version changing at each reboot). A job killed
  this way leaves no result. `scripts/loom_smoke.sh status` reports
  "JOB DIED ... host up since ..." and `launch` replaces the stale marker;
  relaunch. Before blaming your job for a reboot, check whether the kernel
  version changed.
- **The GPU is shared.** Check for other users' GPU processes before timing,
  and wait; never record a baseline under contention.

## Local GPU box (fresh Ubuntu VM with ROCm)

- ROCm is usually at `/opt/rocm` (`ROCM_PATH`), with `hipcc` in
  `/opt/rocm/bin` and `clang`, `llvm-readelf`, `amdgpu-arch` in
  `/opt/rocm/llvm/bin`; neither is on `PATH` by default. The scripts default to
  `${ROCM_PATH:-/opt/rocm}` and detect the GPU target with `amdgpu-arch`
  (fallback `rocminfo`); use that one target for `AMD_ARCH` and every
  `loom-compile --target`, or `hipModuleLoad` fails.
- Your user must be in the `render` and `video` groups to open the GPU
  (`rocminfo` errors otherwise).
- Work dirs default to `$LOOM_WORK`, else `$HOME/loom-work` (`LOOM_REMOTE_BASE` overrides).
- **Every clang uses the newest GCC installed.** After Setup's `g++-14`,
  ROCm's clang (hipcc, and the Loom tool build) uses libstdc++ 14 headers.
  Override per build with `--gcc-install-dir=/usr/lib/gcc/x86_64-linux-gnu/<N>`
  (for hipcc: `HIPCC_COMPILE_FLAGS_APPEND` and `HIPCC_LINK_FLAGS_APPEND`).
  Ubuntu 24.04 with ROCm already has libstdc++ 11 headers (a ROCm dependency).
- **std::array<T, 0> breaks stock HIP builds (Microphysics #2115).** For a
  network without intermediate rates (`triple_alpha_plus_cago`,
  `ignition_simple`, `powerlaw`, ...), `locate_intermediate_rate_index`
  (`networks/rhs.H`) indexes a `std::array<int, 0>`; with libstdc++ >= 13
  its `operator[]` is a host-only function and hipcc fails:
  `reference to __host__ function 'operator[]' in __host__ __device__ function`.
  Fixed upstream by PR #2116; until that is in your Microphysics, either apply
  it locally or build against libstdc++ 11, which compiles it (GCC 13 doesn't):
  `export HIPCC_COMPILE_FLAGS_APPEND=--gcc-install-dir=/usr/lib/gcc/x86_64-linux-gnu/11
  HIPCC_LINK_FLAGS_APPEND=--gcc-install-dir=/usr/lib/gcc/x86_64-linux-gnu/11`
  before `ship_test_react.sh --build`. Record which one in the agent note.

## Remote host facts (moth, gfx90a MI210)

Only when the GPU is on a remote host: set `LOOM_HOST=moth
LOOM_REMOTE_BASE=/mnt/data/bwibking` (and `ROCM_PATH=/opt/rocm/core-10.1`
for the scripts that run there); the defaults are for a local GPU box.

- Work under `/mnt/data/bwibking/` (home is a full NFS quota).
- ROCm 10.1 at `/opt/rocm/core-10.1`; set
  `IREE_HAL_AMDGPU_LIBHSA_PATH=/opt/rocm/core-10.1/lib/libhsa-runtime64.so.1`
  for HRX, or it picks up an older libhsa.
- `loom-import-cxx` needs C++23 libstdc++: build with conda GCC 14 at
  `/mnt/data/bwibking/gcc14` (`--gcc-install-dir=...`, rpath to its lib),
  `-DLOOM_IMPORT_CXX=ON`; tools from `/mnt/data/bwibking/venv-py312/bin`.
- Ship source with `git archive` piped over ssh into a staging dir, then
  `rsync -rlpc --delete --exclude=/build/`. Never repack with macOS `tar`
  (AppleDouble `._*` files broke CMake configure); if you must, set
  `COPYFILE_DISABLE=1`.
- Don't multiplex two payloads through one ssh stdin. A trailing `&`
  backgrounds the *whole* `&&` list, and a background job's stdin is
  `/dev/null`: `ssh host 'tar xzf - && cd d && nohup job &' < payload.tgz`
  gives tar an empty stream (it happened again on 2026-10-18). Unpack in the
  foreground, then background only the job on its own line:
  `tar xzf - -C dir || exit 1` newline `nohup setsid bash job.sh > log 2>&1 < /dev/null &`.
- Small payloads can go over ssh stdin as a tarball; build it with
  `COPYFILE_DISABLE=1 tar --no-xattrs czf ...` on macOS, or GNU tar on the
  remote end prints a warning per file for `com.apple.provenance` headers.

## Shell and scripts

- **Use bash scripts for anything with flags in variables.** The user's shell
  is zsh, which does not word-split `$FLAGS` and expands `=word`; this broke
  importer probes and an orchestrator's cutoff check. Put waits, retries and
  chained remote steps in a `bash` script file, not inline.
- **`pgrep -f` / `pkill -f` match their own command line.** A wait loop
  `while pgrep -f loomval_burn; do sleep 30; done` run through `bash -lc
  '...loomval_burn...'` finds itself and never ends (cost 8 min in the E2E
  run). The `[l]oomval_burn` bracket trick only helps when the name appears
  nowhere else in the issuing command: `until pgrep -f "[l]oomval_sens" &&
  test -e .../loomval_sens3d.gnu.ex` still matched itself through the path
  (10 more min in the second run). Wait on a PID instead (`wait $pid`,
  `kill -0 $(cat pidfile)`), or use `loom_smoke.sh wait` for Loom compiles.
- **Never estimate the time.** An agent guessed 13:25 when it was 12:07 and
  quit 2 h early. Run `date` before any time-based decision or timestamp.
- After an interrupted command, check for partial state (a rejected
  init+setup still created a worktree and state dir once).

## Correctness traps

- **Validate where the hard paths are.** The 5-step check never reached the
  perturbation step, so for 9 campaigns perturbed inputs were never tested. For test_react
  that means the full T/ρ/X grid (`n_cell` ≥ 16), including zones that
  fail, retry, hit NSE, or take many steps. Compare step counts too.
- **GCC constant-folds libm calls in validators.** With `-O3`, g++ peeled
  the first iteration of a validator's `for (lt = 4.0; ...) T = pow(10, lt)`
  loop and evaluated the stock side's `cbrt(1e-4)` at compile time with MPFR
  (correctly rounded), while the Loom side called glibc's `cbrt` at run time
  (1 ulp off): 2,104 false 1-ulp mismatches, all at T = 1e4. Keep each
  per-sample validator function `__attribute__((noinline))` (the shipped
  templates do) or feed inputs through `volatile`. A standalone test that
  gives the same value on both sides is the tell.
- **Constant-argument libm calls on the GPU.** hipcc can fold `exp(24.1)` at
  compile time while Loom computes it at run time; see
  `microphysics-translation.md`, "Math", for the Phase 6 check.
- **FP contraction in reference builds.** hipcc
  (`-ffp-contract=fast-honor-pragmas`) and GCC fuse `a*b + c` by default;
  policy and flag in SKILL.md, "Fixed policies".
- **Decide the comparison target before looking at results** (bit-exact
  here), and never relax it to pass. Non-bit-exact rewrites (reassociated sums, a forward-differentiated
  Jacobian) drifted trace abundances by ~1e-14 over 50 steps; the user
  approved each one and the tolerance change explicitly.
- **Kernel identity and layout.** The harness refuses a code object whose
  entry points, parameters, workgroup size, or record layout differ from the
  host build's. Do the same: `static_assert` sizeof/offsetof on both sides.
- **Workgroup count is a compile-time config.** `loom-compile` binds
  `<root>.workgroup_count.x` at compile time. Changing `n_cell` or
  `max_grid_size` without recompiling launches the wrong grid. Guard
  `idx < ncells` in the kernel and have the host check the count it was
  compiled for (pass it in the params, compare at launch).
- **Runtime choice.** Record which `libamdhip64` (ROCm or HRX) each result
  used. HRX's runtime does not support managed memory, and test_react's
  runtime parameters are `AMREX_GPU_MANAGED`, so the AMReX side of
  test_react should run on ROCm's runtime unless you have removed managed
  memory.
- **hipcc code objects.** `hipcc --genco` emits an offload bundle; for a bare
  ELF use `--cuda-device-only --no-gpu-bundle-output`.

## Frontend and compiler traps (re-probe; some may be fixed)

- `return` inside a loop fails (LOWERING/059). Use a result code and `break`.
- SUBRANGE/026 at compile time means a missing bounds assumption:
  `MPL_ASSUME` (`loom-frontend-rules.md`, SUBRANGE row). If the assumption
  is there and the bound in the message is the whole `short` range
  ("may reach 262152"), the tools predate hrx-system #2135 (issue #2122,
  probe f20).
- A `return` inside an `if` nested in another `if`, with code after both
  (stock `actual_log_screen`'s Debye-Hückel skip), fails import:
  "conditional returns require at most one fallthrough arm" (probe f24; a
  single-level early return imports, f10, f25). Keep a `done` flag there.
- Uninitialized struct locals (`burn_t s;`) fail; write `T x{}`. Scalar and
  array locals may stay uninitialized (hrx-system #2105).
- Text Loom modules over 65,535 lines were rejected (#1877).
- BACKEND/005 `spill-traffic-register-exhausted` can appear from harmless-looking
  statement moves (#1879, #2016). Keep the failing source, reduce it, report it.
- `rocprofv3` misreports VGPR/AGPR columns on gfx90a (ROCm/rocm-systems#13686).
  Read `.vgpr_count`/`.agpr_count` from `llvm-readelf --notes` or the kernel
  descriptor instead.
- AGPRs share the 512-entry file with VGPRs on gfx90a/gfx942; "spill to AGPR"
  adds no capacity at 2 waves/SIMD.

- **Host validation can't see stock names.** Validators include the stock
  headers, so a translated function that still calls `autodiff::eval`,
  `amrex::...`, a `*_rp::` parameter or another stock name compiles and
  validates bit for bit on the host (it happened in `rate_al27pg`). The
  `no-stock-namespace` lint rule and the Loom smoke import catch it.
- **A full fused burn takes ~30 min to compile, almost all in register
  repair.** There is no `noinline` from C++, so EOS + rates + screening +
  Jacobian + integrator land in one function (aprox13 + helmholtz + VODE:
  924k values growing to 1.33M). The spill-free frame build (`frame.c`) then
  loops schedule → allocate → rematerialize/split/materialize spills, each
  round a whole-function schedule + allocation; the round budget is the value
  count, so it never times out on its own. Small module smokes take under 1 min.
- **Watch the repair loop, don't guess.** A long compile has no output until it
  finishes. The compiler-diagnostic tools in the `loom-burn-kernel-tools` follow-up print one line per
  register-repair round and sample the compiler's stack.
- **Check scratch against the hardware limit before any GPU run.** gfx90a
  allows at most 131,056 B of private segment per lane ((2^13 - 1) KiB per
  wave64). hipcc refuses to link anything larger. `loom-compile` does too (hrx-system
  #2091), at the very end of the compile: AMDGPU/053 "requires N private bytes
  per workitem, but processor 'gfx90a' supports at most 131056". (Tools
  without #2091 write the code object anyway and the launch fails with
  `HSA_STATUS_ERROR_OUT_OF_RESOURCES`.) aprox13 + helmholtz + VODE fits at the
  pin (about 72 KB).
- **A smoke root's code object can't stand in for the production root.** The
  arguments, per-zone I/O and baked `workgroup_count` differ, so test_react
  can't launch it. Once the burn body compiles, start the production root
  compile (`loom_gpu_compare.sh`) right away, in parallel with any other
  investigation, rather than paying the ~30 min twice in series.
- **Smoke-import each module as soon as it validates**, not at the end: the
  first wave-1 import found a construct no probe covered (struct of table
  pointers read through a pointer), and the network module compiled only
  after adding `MPL_ASSUME`s. Use `scripts/loom_smoke.sh`.

## Bookkeeping traps

- Record the hrx-system commit and tool md5 with every result. Compiler
  merges moved this kernel by −13% and +10% with no source change; results
  across commits are not comparable.
- Keep every candidate's exact source (commit or patch) next to its result.
  A result with no reproducible source is worthless later.
- The `/tmp` scratchpad is wiped on reboot. Copy validators and reducers you
  want to keep into the agent note's directory.
- When delta-debugging a compiler failure, drop everything outside the keep
  set each round (a reducer that re-added earlier deletions gave a false
  minimal case).
