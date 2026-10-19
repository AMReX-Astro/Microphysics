# Performance work (only after GPU correctness passes)

From 12 campaigns on the chemistry advance kernel (MI210, gfx90a): Loom went
from 2.41× slower than the original HIP kernel to about 0.86–0.89× of it,
2.7× faster than the seed, by editing only the kernel source. These are
hypotheses for a new kernel, not rules; the compiler's response to source
edits was often counter-intuitive.

## Measurement rules

- **Freeze the compiler.** Build the tools once from one hrx-system commit,
  check their md5 on every run. Compiler merges alone moved this kernel by
  −13% and +10%.
- **Paired, interleaved timing on one GPU.** Alternate candidate and best
  (A B B A ...), at least 5 pairs; speedup per pair = best / candidate.
  Keep a change only if the median beats the bar and every clean pair agrees
  in direction (one dissent allowed at ≥ 4 pairs). Set the bar from an A/A
  baseline: max(1.02, 1 + 2 × median A/A deviation). Never set it from the
  maximum spread.
- **Time the burn kernel with events**, after a warm-up, not test_react's
  wall-clock "Run time" (it includes host work and the plotfile).
- **Telemetry with every timed run** (`amd-smi`: GFX clock, power, HBM
  temperature, utilization, other processes). Re-run pairs whose busy clock
  falls below 97% of the median. HIP's chemistry kernel was bimodal
  (~183.5 vs ~191.5 ms) and the slow state tracked HBM temperature ≥ 86 °C
  (DRAM refresh doubling, inferred); Loom's time didn't. Report the HIP state
  when quoting a Loom/HIP ratio.
- **Pin the host side**: timing process on the GPU's NUMA node, builds on
  other cores. Pin GPUs with `ROCR_VISIBLE_DEVICES` only
  (`HIP_VISIBLE_DEVICES` on top masks twice).
- **Matched-source comparison.** A Loom/HIP ratio against the *original*
  Microphysics code compares two sources as much as two compilers. Also
  compile the Loom kernel source with hipcc (an adapter header for the
  `__loom__` bindings) with both OCML and the same f64 recipes. In the
  campaigns this changed the story twice: 4.2% behind hipcc at one point,
  12.9% ahead later.

## Static screens

VGPR/SGPR/scratch/instruction counts predict timing only loosely (rank
correlation about −0.55 over 111 candidates; 5 of 19 kept wins had *worse*
stats). Use a screen to skip only compile failures and byte-identical code
objects. Always time loop or control-structure changes (LU, solve, stage
loop, unrolling, register vs private-array placement). Drop a straight-line
reorder only if it is far worse: scratch > +1.5 KB/lane, VGPRs > 256, or
instructions > +5%. Write down why each candidate was dropped.

## Where the time went

- **Private-memory round trips dominate.** The compiler drains `vmcnt(0)`
  before a reload of private memory, so every spill reload and every load
  after a store to the same private array stalls. Most wins removed a
  store→reload pair from hot code.
- **Register pressure caps occupancy.** The kernel sat at 256 VGPRs (2
  waves/SIMD on gfx90a). Builds that dropped to 248–254 VGPRs while adding
  scratch traffic got *slower*. AGPRs add no capacity at 2 waves.
- **Hoisted constants.** ~124 f64 rate constants hoisted by LICM were
  materialized and spilled in the preheader; ~37 VGPRs held exp/log recipe
  constants for the whole kernel.
- **The Jacobian was ~63% of time; LU + solve ~30% of PC samples.**

## Levers that won (chemistry kernel, roughly in impact order)

1. Load-all / compute / store-all instead of load-update-store loops on
   private arrays (up to +15.8% for one LU change).
2. Fused/pipelined LU with a batched solve; LAPACK-form full-row pivoting so
   the solve uses constant indices held in registers (+6.2%).
3. Forward-differentiated Jacobian blocks with SSA accumulators (+6.9%;
   *not* bit-exact, user-approved).
4. Counters (`n_jac`, `n_step`, `n_accept`) derived or stored after the loop
   instead of updated in hot code (+5.5% for one counter).
5. Tolerances and controller settings read once into constants for the
   error norm (part of +9.3%).
6. Order RHS sub-computations to minimize live values across calls (+2.7%).
7. Accumulate large output sums into running SSA partials as terms become
   ready (adv-mi210-8: two kept runs, +5.5% and +3.0%, ~1.10× combined;
   reassociation, not bit-exact, user-approved).
8. 32-bit `int` loop indices with compile-time trip counts instead of 64-bit
   `size_type`.

## FMA contraction (optional, this phase only)

Every correctness reference is built with `-ffp-contract=off`. Here, and
only here, FMA may be tried as its own candidate: enable it on both sides
being compared (hipcc `-ffp-contract=fast` or explicit `fma` calls in the
kernel source; whatever Loom offers on its side), time it paired like any
other change, and check it against the contraction-off references with
tolerances the user approves for it. Record it as a non-bit-exact change.

## Levers that lost

- Fully unrolled LU (−7 to −13% every time); stage-loop unroll (−24%).
- `+=` on mutable scalars instead of `const` chains (−12%).
- Moving integrator state between registers and memory (−9 to −25%).
- Pushing exp/log constants out of registers (−10 to −40%).
- Select-form choice helpers, sinking arm-only work into branches, merging
  same-condition calls (−2.6 to −4.2%).
- Reassociating *intermediate* definitions (−10%), unlike output sums.
- Stacking several +0.4% wins on top of a +2.7% one (net −2.7%): effects are
  strongly non-additive; re-time stacks, don't sum them.
- Edits to cold code (input checks) reshuffled allocation by ±10–25%.

## Process

- One commit per experiment, one experiment at a time, the exact source of
  each result kept (tag discarded commits).
- Time-box each idea to two attempts or ~30 min and write a verdict.
- Re-profile after each new best and name the hot window an idea targets.
- If running an unattended loop, the campaigns' hrx-system `autoresearch/`
  harness (`ar`, `remote.sh`, `program.md`) is not public (it lived on the
  author's `react-kernels` branch); a new loop should copy its rules: its rules (edit only the kernel dir, frozen tools, interface check,
  correctness before timing, clock line, no local Loom compiles) encode
  most of the footguns. Only start one with the user's go-ahead.
