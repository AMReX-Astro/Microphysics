# Local hrx-system patches

Local patches, not upstream. Both apply cleanly to public hrx-system main
`31bf125137` (the current pin, checked 2026-10-20) and to the previous pin
`4031286f0b` (= fork commit `7cf21edd85`, where they were written). Apply at the repository root with `patch -p1 < FILE`. The usual way is to let
`scripts/build_loom_tools.sh` copy the source tree and apply them:

```bash
SRC=/path/to/hrx-system OUT=/path/to/tools-build \
PATCHES="/path/to/patches/alloc-relocation-overlap.patch /path/to/patches/trace-repair.patch" \
nohup setsid bash build_loom_tools.sh > build.log 2>&1 < /dev/null &
```

## alloc-relocation-overlap.patch

`loom/src/loom/codegen/low/allocation/interval_assignment.c`,
`try_relocate_independent_assignment`: precomputes the lifetime-overlapping
blockers once per relocation into a bitmap over candidate bases (up to 1024
bases; above that it falls back to the original per-base loop over all
assignments). It visits the same bases in the same order with the same
checks, so it picks the same base: code objects are byte-identical with and
without the patch (verified on the smoke eos, math, rates and network
kernels with `scripts/compare_loom_tools.sh`).

Why: each relocation was O(bases x assignments). The full fused burn kernel
was still in allocation after more than 27 minutes without the patch; with it,
the whole compile takes about 32 minutes.

## trace-repair.patch

`loom/src/loom/codegen/low/frame.c`: with `LOOM_TRACE_REPAIR=1` in the
environment, `loom-compile` prints one line per spill-free repair round on
stderr:

```
LOOM_REPAIR round=N limit=N values=N schedule_errors=N allocation_errors=N spill_plans=N rematerialized_operands=N split_operands=N round_seconds=S
```

`round=0` begins a new spill batch; `spill_plans=0` means the repair
converged. Off by default, and output is unchanged without it. Summarize a
trace with `scripts/repair_summary.py`.

Note: hrx-system PR #2047 (scheduler changes) did not change these repair
counters for the fused burn kernel.
