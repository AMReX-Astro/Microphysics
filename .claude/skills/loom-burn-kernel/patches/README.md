# Local hrx-system patches

Local patches, not upstream, against public hrx-system main `fada6788a1`
(the current pin, checked 2026-10-20). Apply at the repository root with
`patch -p1 < FILE`. The usual way is to let `scripts/build_loom_tools.sh`
copy the source tree and apply them:

```bash
SRC=/path/to/hrx-system OUT=/path/to/tools-build \
PATCHES="/path/to/patches/assume-cast-chain.patch /path/to/patches/trace-repair.patch" \
nohup setsid bash build_loom_tools.sh > build.log 2>&1 < /dev/null &
```

`alloc-relocation-overlap.patch` (needed up to `31bf125137` to finish a
full-burn compile) is gone: hrx-system #2107 (`b0ae10b7b6`) filters the
assignment history once per relocation upstream, the same fix. The aprox13
full burn compiles without it (see `assets/README.md`).

## assume-cast-chain.patch (required)

`loom/src/loom/transforms/{scalar,vector}/combine_patterns.c`: compose a
cast chain only with the conversion that directly defines the outer cast's
input. Since hrx-system #2074 (`68905de218`, bisected 2026-10-20) the combine
found the inner cast through `loom_rewriter_value_defining_op`, which looks
through value identities such as `scalar.assume`. For an index loaded from a
`short` field, `extsi i16->i32` → `scalar.assume [1, 6]` → (accessor)
`extsi i32->i64` then folded into one `extsi i16->i64` and the assumed range
was lost, so `loom-compile` failed SUBRANGE/026 ("access byte end may reach
262152", i.e. the whole `short` range) on VODE's `el(L)`, `tau(NQ)`. The
imported IR is identical before and after #2074; the bug is in compile.
Reproducer: probe `f20_assume_short_field_member_index` (fails without the
patch, compiles with it). Probe results are otherwise identical with and
without the patch. Reported upstream as hrx-system #2122; drop the patch
once a fix is in the pinned commit (rerun f20 without it).

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
