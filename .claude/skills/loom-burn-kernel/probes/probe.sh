#!/bin/bash
# Phase 1 Loom frontend probes for a Microphysics burn port (loom-burn-kernel skill).
# One construct per file in src/ (root probe::k), compiled with ONE frozen tool set
# built from the pinned hrx-system commit. Writes probe-summary.txt (tool md5, hrx
# commit, per-probe result, VGPR/AGPR/scratch, and for m_* probes the f64
# instructions emitted) for pasting into the port's agent note.
#
# Usage: probe.sh          (run on the GPU host from anywhere; configuration via environment)
#        probe.sh -h       (this text)
#
# Copy this directory into the port's agent-notes dir, add probes for anything
# the inventory found that src/ doesn't cover, put the skill's OCML recipe
# header next to it as include/f64_math.h (the m_f64_math_h_recipes probe
# includes "f64_math.h"):
#   mkdir -p <probes dir>/include
#   cp <skill>/assets/kernel/f64_math.H <probes dir>/include/f64_math.h
# then run it on the GPU host:
#   W=<probes dir> TOOLS=... bash <probes dir>/probe.sh > probe-run.log 2>&1
#
# Check m_fp_contract: it must show v_mul_f64 + v_add_f64 and no v_fma_f64.
#
# Environment:
#   TOOLS             dir with loom-import-cxx and loom-compile built from the
#                     pinned hrx-system commit                        (required)
#   LOOM_WORK         base dir for the defaults below                 (default $HOME/loom-work)
#   W                 the probes dir (contains src/; writable)        (default: this script's dir)
#   F64_INC           dir with f64_math.h for m_f64_math_h_recipes    (default $W/include; falls
#                     back to ../assets/kernel/f64_math.H if the whole skill was shipped)
#   HRX               hrx-system checkout, only to record its commit  (default $LOOM_WORK/hrx-system;
#                     $TOOLS/SOURCE from build_loom_tools.sh is used first)
#   TARGET            amdgpu target     (default: local GPU from amdgpu-arch/rocminfo, else gfx90a)
#   CPUS              taskset CPU list away from timing cores         (default: all CPUs)
#   ROCM_PATH         ROCm install                                    (default /opt/rocm)
#   L                 dir with llvm-readelf, llvm-objdump (default $ROCM_PATH/llvm/bin, else
#                     $ROCM_PATH/lib/llvm/bin)
#   TMPDIR            kept if set, else $LOOM_WORK/tmp
set -u
usage() { sed -n '2,/^set -u/p' "$0" | sed '$d'; exit "${1:-2}"; }
case "${1:-}" in -h|--help) usage 0;; '') ;; *) usage;; esac
[ -n "${TOOLS:-}" ] || { echo "set TOOLS to a dir with loom-import-cxx and loom-compile built from the pinned commit"; usage; }
here=$(cd "$(dirname "$0")" && pwd)
R=${LOOM_WORK:-$HOME/loom-work}
# default: the kit this script lives in (the copy in the port's agent-notes
# dir, SKILL.md Phase 1); it must be writable (out/, probe-summary.txt)
W=${W:-$here}
HRX=${HRX:-$R/hrx-system}
ROCM_PATH=${ROCM_PATH:-/opt/rocm}
L=${L:-}
if [ -z "$L" ]; then L=$ROCM_PATH/llvm/bin; [ -x "$L/llvm-readelf" ] || L=$ROCM_PATH/lib/llvm/bin; fi
[ -x "$L/llvm-readelf" ] || { echo "no llvm-readelf in $L (set L or ROCM_PATH)"; exit 2; }
if [ -z "${TARGET:-}" ]; then
    TARGET=$({ "$ROCM_PATH/llvm/bin/amdgpu-arch" || "$ROCM_PATH/bin/rocminfo"; } 2>/dev/null | grep -m1 -o 'gfx[0-9a-f]*')
fi
[ -n "$TARGET" ] || { echo "WARNING: no TARGET given and no GPU detected; falling back to gfx90a" >&2; TARGET=gfx90a; }
CPUS=${CPUS:-0-$(( $(nproc) - 1 ))}
export TMPDIR=${TMPDIR:-$R/tmp}
mkdir -p "$TMPDIR"
cd "$W" || exit 1
F64_INC=${F64_INC:-$W/include}
if [ ! -f "$F64_INC/f64_math.h" ] && [ -f "$here/../assets/kernel/f64_math.H" ]; then
  # whole skill shipped: stage its recipe header under the name the probe includes
  F64_INC=$W/out/f64_inc
  mkdir -p "$F64_INC" && cp "$here/../assets/kernel/f64_math.H" "$F64_INC/f64_math.h"
fi

for t in loom-import-cxx loom-compile; do
  [[ -x $TOOLS/$t ]] || { echo "MISSING $TOOLS/$t"; exit 2; }
done

OUT=out
mkdir -p "$OUT"
SUMMARY=probe-summary.txt
{
  echo "# date: $(date)"
  echo "# host: $(hostname)"
  echo "# hrx-system: $(sed -n 's/^commit: //p' "$TOOLS/SOURCE" 2>/dev/null | grep . || git -C "$HRX" rev-parse HEAD 2>/dev/null || cat "$HRX/COMMIT" 2>/dev/null || echo unknown)"
  echo "# f64_math.h: $([ -f "$F64_INC/f64_math.h" ] && md5sum "$F64_INC/f64_math.h" || echo "MISSING in $F64_INC (m_f64_math_h_recipes will fail)")"
  echo "# tools md5:"
  md5sum "$TOOLS/loom-import-cxx" "$TOOLS/loom-compile" | sed 's/^/#   /'
  echo "# target: $TARGET"
  printf '%-34s %-13s %s\n' probe result detail
} > "$SUMMARY"

for f in src/*.cpp; do
  p=$(basename "$f" .cpp)
  o=$OUT/$p
  extra=()
  [[ $p == m_f64_math_h_recipes ]] && extra=(--I="$F64_INC")
  if ! nice -n 10 taskset -c "$CPUS" timeout 300 "$TOOLS/loom-import-cxx" --data-model=lp64 --fp-model=precise \
        --root=probe::k "${extra[@]}" --output="$o.loom" "$f" > "$o.import.log" 2>&1; then
    printf '%-34s %-13s %s\n' "$p" import-fail \
      "$(grep -m1 -E 'error|ERROR|[A-Z]+/[0-9]{3}' "$o.import.log" | sed 's/^[^ ]*: //' | cut -c1-160)" >> "$SUMMARY"
    continue
  fi
  if ! nice -n 10 taskset -c "$CPUS" timeout 600 "$TOOLS/loom-compile" "$o.loom" --root=probe.k \
        --config=probe.k.workgroup_count.x=1 --config=probe.k.workgroup_count.y=1 --config=probe.k.workgroup_count.z=1 \
        --format=amdgpu-hsaco --target=amdgpu:"$TARGET" --compile-report=summary \
        --compile-report-output="$o.report.json" --output="$o.hsaco" > "$o.compile.log" 2>&1; then
    printf '%-34s %-13s %s\n' "$p" compile-fail \
      "$(grep -m1 -E 'error|ERROR|RESOURCE|[A-Z]+/[0-9]{3}' "$o.compile.log" | cut -c1-160)" >> "$SUMMARY"
    continue
  fi
  n=$("$L/llvm-readelf" --notes "$o.hsaco" | awk '{for (j = 1; j < NF; j++) { if ($j == ".vgpr_count:") v = $(j+1); if ($j == ".agpr_count:") a = $(j+1); if ($j == ".private_segment_fixed_size:") s = $(j+1) }} END{printf "vgpr=%s agpr=%s scratch=%sB", v, a, s}')
  # For math probes, record which f64 instructions the op lowered to (native vs. expanded).
  if [[ $p == m_* ]]; then
    n="$n | $("$L/llvm-objdump" -d "$o.hsaco" | grep -oE '^\s+v_[a-z0-9_]+(f64|f32)' | sort | uniq -c | sort -rn | head -6 | awk '{printf "%s:%s ", $2, $1}')"
  fi
  printf '%-34s %-13s %s\n' "$p" ok "$n" >> "$SUMMARY"
done
echo "# done: $(date)" >> "$SUMMARY"
cat "$SUMMARY"
