# Hopper Primordial Rodas5P Prototype Plan

## Summary

Build a CUDA/Hopper-only prototype for `primordial_chem + Rosenbrock + Rodas5P` that bypasses the current per-zone `burner()` path and launches a cooperative block kernel over batches of zones. The performance hypothesis is narrower than tensor-core/TMA warp-specialization papers such as Tawa: this kernel is CUDA-core-heavy chemistry plus a small dense solve, so the first prototype optimizes for reduced local-memory spills, improved scratch locality, and explicit inter-warp work partitioning, not Tensor Core utilization or immediate production integration.

Default design choices:

- Target NVIDIA Hopper/GH200 only, compiled for `sm_90`.
- Support only Strang, primordial chemistry, analytic Jacobian, number-density mode, and `integrator.rosenbrock_tableau = 0` / Rodas5P.
- Use `16 zones/block`, `8 warps/block`, one cooperative CTA per zone tile.
- Use dynamic shared memory with explicit Hopper opt-in and preferred shared-memory carveout.
- Treat Tawa-style asynchronous dataflow as a design checklist, not as the v1 implementation model. The prototype should first prove correctness and spill reduction with simple barriers, then only add mbarrier/named-barrier channels when profiling shows they are needed.
- Keep all existing generic Rosenbrock/VODE paths unchanged unless the new prototype test explicitly calls the Hopper kernel.

## Tawa-Informed Design Constraints

- Do not assume Tawa's GEMM/attention speedups transfer directly. Tawa primarily overlaps TMA global-to-shared transfers with WGMMA Tensor Core work; this prototype should instead measure whether cooperative CUDA-core execution reduces local-memory traffic enough to matter.
- Preserve a first-class communication model between warp roles. Even if v1 lowers this to `__syncthreads()`, design the code and generator around explicit producer/consumer phases so a later mbarrier/named-barrier implementation is mechanical rather than a rewrite.
- Keep pipeline depth and buffering as measured parameters. Tawa's deeper buffers help until shared memory and registers reduce occupancy; this prototype already starts near a large shared-memory footprint, so `16 zones/block` and scratch size must be validated with ptxas and Nsight data.
- Prefer generated role partitions over hand-written expression movement. The RHS/Jacobian expression graph is too large to maintain safely by hand, and any partitioning must preserve expression order unless a later validation step proves an algebraic rewrite is acceptable.

## Key Implementation Changes

### Prototype Entry Point

- Add a new Hopper-specific batch API, separate from `burner()`:
  - `primordial_rodas5p_hopper_burn(burn_t* states, amrex::Real* dt, int nstates)`
  - CUDA-only host launcher plus device kernel.
- Do not call this from `actual_integrator.H` in v1. The current `actual_integrator()` remains the scalar per-zone reference path.
- Add a dedicated prototype test executable that initializes an array of `burn_t` states, runs both:
  - existing scalar `burner(state, dt)` reference,
  - new Hopper batch kernel,
  then compares results.

### Kernel Structure

- One CUDA block owns a tile of `16` zones.
- Each warp specializes in part of the Rodas5P work:
  - warp 0: state setup, validation, timestep controller, reductions, success/error bookkeeping.
  - warps 1-3: RHS expression chunks.
  - warps 4-6: analytic Jacobian expression chunks and Jacobian row/block assembly.
  - warp 7: shared-memory LU/solve coordination, with helpers called by all warps where useful.
- Use `__syncthreads()` for v1 synchronization. Avoid Hopper named barriers and thread-block clusters in the first prototype, but keep phase boundaries explicit enough to replace whole-CTA barriers with targeted producer/consumer handshakes later.
- Hard-code Rodas5P coefficients and unroll all 8 stages in the prototype kernel. Do not use runtime tableau dispatch.

Synchronization milestones:

- Milestone 1: whole-CTA barriers only, with one shared-memory layout and deterministic compare against the scalar reference.
- Milestone 2: if Nsight shows barrier stalls or idle specialized warps dominate, introduce aref-like double-buffered channels for RHS/Jacobian/LU handoff using named barriers or mbarriers.
- Milestone 3: only after Milestone 2 data, consider Hopper-specific features such as thread-block clusters or TMA-style movement for any bulk global/shared transfers that remain on the critical path.

### Shared Memory Layout

Use one dynamic shared-memory allocation per block with explicit byte offsets and padding/alignment:

- `y[15][16]`, `ynew[15][16]`, `work[15][16]`.
- `ak[8][15][16]` for Rodas5P stages.
- `jac[15][15][16]`, factored in place as `A = 1/(h*gamma) I - J`.
- `rhs_tmp[15][16]`.
- `x_scratch[534][16]` for primordial generated expression temporaries during RHS/Jacobian phases.
- `pivot[15][16]`, status flags, step counters, and small reduction buffers.

Expected shared-memory budget:

- `x_scratch`: about `68.4 KiB`.
- `jac`: about `28.8 KiB`.
- Rodas5P stage/state/RHS buffers: about `25-30 KiB`.
- Total target: `<= 140 KiB` before padding, comfortably below Hopper's 227 KiB per-block limit.

This budget should be treated as a first correctness milestone, not the final performance layout. The generator should eventually compute liveness ranges for `x_scratch` temporaries and reuse slots where possible. Before trying `32 zones/block` or deeper buffering, measure active CTAs/SM, register pressure, spill load/store counts, and shared-memory bank conflicts with the uncompressed layout.

### Generated Primordial RHS/Jacobian Prototype

- Add a small generator script for the prototype instead of hand-maintaining thousands of expressions.
- Input: existing `networks/primordial_chem/actual_rhs.H`.
- Parse top-level generated assignments in:
  - `rhs_specie(...)`,
  - `rhs_eint(...)`,
  - `jac_nuc(...)`.
- Emit a Hopper-only generated header with:
  - phase functions that write `xN` temporaries into `x_scratch`,
  - RHS assembly into shared `rhs_tmp`,
  - Jacobian assembly into shared `jac`.
- Preserve expression order exactly for v1. Partition only by contiguous statement ranges and Jacobian output blocks; do not algebraically rewrite expressions.
- If parsing fails on an unexpected statement form, fail code generation loudly rather than silently changing math.
- After correctness is established, add a generator mode that performs liveness analysis for generated temporaries and reports the peak scratch footprint before emitting any compacted shared-memory layout.

### Hopper Launch Configuration

- Launcher sets:
  - `cudaFuncAttributeMaxDynamicSharedMemorySize` to the computed shared-memory size.
  - `cudaFuncAttributePreferredSharedMemoryCarveout` to maximum/shared-preferred.
- Kernel launch:
  - `blockDim = 8 * 32 = 256 threads`.
  - `gridDim = ceil(nstates / 16)`.
  - dynamic shared-memory size computed by a constexpr layout helper.
- Compile only when `AMREX_USE_CUDA` is enabled. Non-CUDA builds should not see the Hopper prototype symbols unless explicitly guarded.

## Test Plan

### Correctness Tests

- Add a prototype unit test with three modes:
  - `reference`: existing scalar `burner()` only.
  - `hopper`: Hopper batch kernel only.
  - `compare`: run both on identical states and compare.
- Test state set:
  - the default `unit_test/burn_cell_primordial_chem/inputs_primordial_chem` state,
  - a small grid of temperature/density/redshift-relevant states,
  - edge states with tiny species number densities near `SMALL_X_SAFE`,
  - multi-step states that force at least one Rodas5P rejection.
- Compare:
  - `success`,
  - `error_code`,
  - `n_rhs`, `n_jac`, `n_step`,
  - final `xn[0:NumSpec]`,
  - final `T`, `e`, and `time`.
- Numerical acceptance:
  - exact match for `success` and `error_code`,
  - `n_step`, `n_rhs`, `n_jac` must match in initial prototype,
  - species/energy relative difference `<= 1e-10` for normal values, absolute difference `<= 1e-30` for tiny number densities.

### Build and Spill Tests

- Build CUDA/Hopper target with verbose ptxas output:
  - CMake CUDA build with `AMReX_GPU_BACKEND=CUDA`,
  - architecture `sm_90`.
- Capture and save:
  - registers per thread,
  - spill stores/loads,
  - shared-memory bytes,
  - local-memory bytes.
  - active CTAs per SM / achieved occupancy,
  - shared-memory bank conflicts,
  - barrier stall cycles,
  - local-memory load/store transactions.
- Acceptance for prototype:
  - Hopper kernel compiles without static shared-memory overflow,
  - dynamic shared-memory opt-in succeeds,
  - ptxas spill count is substantially below the existing per-zone kernel baseline,
  - occupancy loss from shared memory and registers is understood and recorded,
  - whole-CTA barrier cost is measured before adding finer-grained synchronization,
  - no unsupported-path fallback is used in the Hopper test.

### Performance Tests

- Benchmark three cases on GH200:
  - existing scalar GPU path if available,
  - existing CPU scalar path as sanity reference,
  - new Hopper batch kernel.
- Measure:
  - wall time for fixed number of zones and timesteps,
  - local load/store transactions,
  - shared-memory throughput,
  - achieved occupancy,
  - branch/replay stalls if available through Nsight Compute.
- Initial performance goal:
  - demonstrate lower local-memory traffic and no correctness regression.
  - wall-time speedup is expected but not required for the first correctness milestone.
- Follow-up performance goal:
  - if the cooperative CTA path is faster but launch or tail effects dominate, prototype a persistent-kernel/work-queue variant with one resident CTA stream per SM and compare it against the simple grid launch.

## Assumptions

- The first prototype is a standalone Hopper batch path, not a production replacement for `burner()`.
- The first target is CUDA/Hopper only; `gfx90a` is intentionally deferred.
- `16 zones/block` is the default because storing all primordial expression scratch plus Jacobian and Rodas5P stage data for `32` zones would leave too little shared-memory headroom.
- The first implementation prioritizes correctness and spill reduction using whole-CTA barriers; named barriers, mbarriers, TMA, clusters, persistent kernels, liveness-compressed scratch, and a `32 zones/block` Hopper-tuned variant are follow-up optimizations gated by measurements.
- Unsupported runtime configurations fall back in tests to the existing scalar path rather than partially running the Hopper prototype.
