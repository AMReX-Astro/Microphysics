#!/bin/bash
# A/B two Loom tool builds (e.g. pinned vs patched, see ../patches/) on the
# same kernel sources: import + compile every spec with both tool sets and
# report compile time, resources and whether the code objects are identical.
#
# Runs on the GPU host (locally, or there over one ssh connection), detached;
# background only the job:
#   TOOLS_A=... TOOLS_B=... KERNEL_DIR=... \
#   nohup setsid bash compare_loom_tools.sh "eos smoke_wave1.cpp mpl::smoke_eos" \
#       "net smoke_burn.cpp mpl::smoke_burn" > ab.log 2>&1 < /dev/null &
#
# Usage: compare_loom_tools.sh [SPEC ...]      (or SPECS env, ';'-separated)
#        compare_loom_tools.sh -h              (this text)
#   SPEC = "NAME SOURCE ROOT": NAME labels the output files, SOURCE is a .cpp
#   relative to KERNEL_DIR, ROOT the C++ root (Loom symbol: ROOT with :: -> .).
#   Each root is compiled with workgroup_count 1x1x1, --fp-model=precise,
#   --data-model=lp64, --compile-report=summary.
#
# Output (in OUT): per spec and tool set NAME.{A,B}.{loom,hsaco,import.log,
# compile.log}, and summary.txt with one line per spec and set:
#   NAME set=A compile=SECONDSs vgpr= agpr= scratch=B md5=...
# followed by "NAME code objects IDENTICAL" or "NAME code objects DIFFER".
# Exit status 0 if every spec compiled with both sets and all are identical.
#
# Environment:
#   TOOLS_A, TOOLS_B  dirs with loom-import-cxx, loom-compile (+ tools.md5)   (required)
#   KERNEL_DIR        dir with the kernel sources (headers next to them)      (required)
#   SPECS             ';'-separated specs, used when no SPEC arguments are given
#   LOOM_REMOTE_BASE  base dir for defaults                     (default $LOOM_WORK, else $HOME/loom-work)
#   OUT               output dir                                (default $LOOM_REMOTE_BASE/loom-ab/<UTC timestamp>)
#   TARGET            amdgpu target       (default: local GPU from amdgpu-arch/rocminfo, else gfx90a)
#   CPUS              taskset CPU list                          (default: all CPUs)
#   ROCM_PATH         ROCm install                              (default /opt/rocm)
#   LLVM              dir with llvm-readelf  (default $ROCM_PATH/llvm/bin, else $ROCM_PATH/lib/llvm/bin)
#   IMPORT_TIMEOUT, COMPILE_TIMEOUT  seconds                    (default 1800, 3600)
#   TMPDIR            temporary dir                             (default $OUT/tmp)
#   Extra environment (e.g. LOOM_TRACE_REPAIR=1) is passed through to the tools.
set -u
usage() { sed -n '2,/^set -u/p' "$0" | sed '$d'; exit "${1:-2}"; }
case "${1:-}" in -h|--help) usage 0;; esac

B=${LOOM_REMOTE_BASE:-${LOOM_WORK:-$HOME/loom-work}}
[ -n "${TOOLS_A:-}" ] && [ -n "${TOOLS_B:-}" ] && [ -n "${KERNEL_DIR:-}" ] || { echo "set TOOLS_A, TOOLS_B and KERNEL_DIR"; usage; }
OUT=${OUT:-$B/loom-ab/$(date -u +%Y%m%dT%H%M%SZ)}
ROCM_PATH=${ROCM_PATH:-/opt/rocm}
LLVM=${LLVM:-}
if [ -z "$LLVM" ]; then LLVM=$ROCM_PATH/llvm/bin; [ -x "$LLVM/llvm-readelf" ] || LLVM=$ROCM_PATH/lib/llvm/bin; fi
[ -x "$LLVM/llvm-readelf" ] || { echo "no llvm-readelf in $LLVM (set LLVM or ROCM_PATH)"; exit 2; }
if [ -z "${TARGET:-}" ]; then
    TARGET=$({ "$ROCM_PATH/llvm/bin/amdgpu-arch" || "$ROCM_PATH/bin/rocminfo"; } 2>/dev/null | grep -m1 -o 'gfx[0-9a-f]*')
fi
[ -n "$TARGET" ] || { echo "WARNING: no TARGET given and no GPU detected; falling back to gfx90a" >&2; TARGET=gfx90a; }
CPUS=${CPUS:-0-$(( $(nproc) - 1 ))}
IMPORT_TIMEOUT=${IMPORT_TIMEOUT:-1800}
COMPILE_TIMEOUT=${COMPILE_TIMEOUT:-3600}

specs=()
if [ $# -gt 0 ]; then
    specs=("$@")
elif [ -n "${SPECS:-}" ]; then
    IFS=';' read -r -a specs <<< "$SPECS"
else
    echo "no specs given"; usage
fi

K=$(cd "$KERNEL_DIR" && pwd) || exit 1
declare -A T=([A]=$TOOLS_A [B]=$TOOLS_B)
for s in A B; do
    for t in loom-import-cxx loom-compile; do
        [ -x "${T[$s]}/$t" ] || { echo "MISSING ${T[$s]}/$t"; exit 2; }
    done
    if [ -f "${T[$s]}/tools.md5" ]; then
        (cd "${T[$s]}" && md5sum -c --quiet tools.md5) || { echo "TOOLS MD5 MISMATCH in ${T[$s]}"; exit 1; }
    fi
done
mkdir -p "$OUT" || exit 1
OUT=$(cd "$OUT" && pwd)
export TMPDIR=${TMPDIR:-$OUT/tmp}
mkdir -p "$TMPDIR"
cd "$OUT" || exit 1
SUM=$OUT/summary.txt
say() { echo "$*" | tee -a "$SUM"; }
: > "$SUM"
say "# $(date) host=$(hostname) target=$TARGET kernel_dir=$K"
for s in A B; do
    say "# $s: ${T[$s]} md5 $(md5sum "${T[$s]}/loom-import-cxx" "${T[$s]}/loom-compile" | awk '{print substr($1, 1, 12)}' | tr '\n' ' ')"
done

# run SET NAME SOURCE ROOT -> sets $md5 (empty on failure)
run() {
    local set=$1 name=$2 src=$3 root=$4
    local o=$name.$set lr=${4//::/.} t0 dt res
    md5=""
    t0=$(date +%s)
    if ! taskset -c "$CPUS" timeout "$IMPORT_TIMEOUT" "${T[$set]}/loom-import-cxx" --data-model=lp64 \
          --fp-model=precise --root="$root" --output="$o.loom" "$K/$src" > "$o.import.log" 2>&1; then
        say "$name set=$set IMPORT FAIL ($(( $(date +%s) - t0 ))s): $(grep -m1 error "$o.import.log" | cut -c1-200)"
        return 1
    fi
    t0=$(date +%s)
    if ! taskset -c "$CPUS" timeout "$COMPILE_TIMEOUT" "${T[$set]}/loom-compile" "$o.loom" --root="$lr" \
          --config="$lr.workgroup_count.x=1" --config="$lr.workgroup_count.y=1" --config="$lr.workgroup_count.z=1" \
          --format=amdgpu-hsaco --target=amdgpu:"$TARGET" --compile-report=summary \
          --compile-report-output="$o.report.json" --output="$o.hsaco" > "$o.compile.log" 2>&1; then
        say "$name set=$set COMPILE FAIL ($(( $(date +%s) - t0 ))s): $(grep -m1 error "$o.compile.log" | cut -c1-200)"
        return 1
    fi
    dt=$(( $(date +%s) - t0 ))
    res=$("$LLVM/llvm-readelf" --notes "$o.hsaco" | awk '{for (j = 1; j < NF; j++) { if ($j == ".vgpr_count:") v = $(j+1); if ($j == ".agpr_count:") a = $(j+1); if ($j == ".private_segment_fixed_size:") s = $(j+1) }} END{printf "vgpr=%s agpr=%s scratch=%sB", v, a, s}')
    md5=$(md5sum "$o.hsaco" | cut -d' ' -f1)
    say "$name set=$set compile=${dt}s $res md5=$md5"
}

status=0
for spec in "${specs[@]}"; do
    read -r name src root extra <<< "$spec"
    if [ -z "${root:-}" ] || [ -n "${extra:-}" ]; then say "bad spec '$spec' (want: NAME SOURCE ROOT)"; status=1; continue; fi
    case "$name" in *[!A-Za-z0-9._-]*) say "bad NAME '$name'"; status=1; continue;; esac
    run A "$name" "$src" "$root"; ma=$md5
    run B "$name" "$src" "$root"; mb=$md5
    if [ -z "$ma" ] || [ -z "$mb" ]; then
        say "$name code objects NOT COMPARED (a compile failed)"; status=1
    elif cmp -s "$name.A.hsaco" "$name.B.hsaco"; then
        say "$name code objects IDENTICAL"
    else
        say "$name code objects DIFFER"; status=1
    fi
done
say "# done $(date)"
exit $status
