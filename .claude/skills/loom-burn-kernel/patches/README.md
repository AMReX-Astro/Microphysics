# Local hrx-system patches

Local patches, not upstream, against public hrx-system main `fa8eba4554`
(the current pin, checked 2026-10-21). Apply at the repository root with
`patch -p1 < FILE`. The usual way is to let `scripts/build_loom_tools.sh`
copy the source tree and apply them:

```bash
SRC=/path/to/hrx-system OUT=/path/to/tools-build \
PATCHES=/path/to/patches/trace-repair.patch \
nohup setsid bash build_loom_tools.sh > build.log 2>&1 < /dev/null &
```

No patch is required at `fa8eba4554`; `trace-repair.patch` is a diagnostic.
It also applies cleanly at main `970d39ac1e` (checked 2026-10-27).
Dropped patches:
- `assume-cast-chain.patch` (required at `fada6788a1`): hrx-system #2135
  (`653fe783ec`) fixes #2122 upstream; producer matching in the combine now
  returns the operand's own SSA definition, so an assumed range on an index
  loaded from a `short` field survives. Probe
  `f20_assume_short_field_member_index` compiles with unpatched tools.
- `alloc-relocation-overlap.patch` (needed up to `31bf125137`): hrx-system
  #2107 (`b0ae10b7b6`) filters the assignment history once per relocation
  upstream, the same fix.

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
