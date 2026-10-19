#!/bin/bash
# Phase 6 on the GPU host: compile the Loom kernel root for the production
# zone count, check its argument metadata, and run test_react in compare mode
# (stock HIP burn vs Loom burn, same process, bit-for-bit) for each run spec.
#
# Run on the GPU host (locally, or there over one ssh connection: unpack in
# the foreground, background only this script; footguns.md), detached:
#   nohup setsid bash loom_gpu_compare.sh > compare.log 2>&1 < /dev/null &
#
# Usage: loom_gpu_compare.sh        (all configuration via environment)
#        loom_gpu_compare.sh -h     (this text; also printed when a required
#                                    variable is missing)
#
# Environment (required):
#   TEST_DIR    test_react directory of the HIP build tree (contains the built *.hip*.ex)
#   KERNEL_SRC  kernel root source, relative to TEST_DIR (e.g. loom/<port>/burn_kernel.cpp)
#   INPUTS      inputs file in TEST_DIR (e.g. inputs_<network>)
#   TOOLS       frozen loom-import-cxx / loom-compile (+ tools.md5) from the pinned commit
# Environment (optional):
#   ROOT        C++ root name (Loom symbol is ROOT with :: -> .)   (default mpl::burn_zones)
#   NCELL       test_react n_cell; zones = NCELL^3, one box (max_grid_size = NCELL) (default 16)
#   RUNS        ';'-separated extra runtime-parameter sets, one compare run each;
#               an empty set = the inputs file alone             (default "": one run)
#               e.g. RUNS=";eos.prad_limiter_rho_c=1.0e7 eos.prad_limiter_delta_rho=3.0e6"
#   TARGET      amdgpu target   (default: local GPU from amdgpu-arch/rocminfo, else gfx90a);
#               must match the AMD_ARCH test_react was built with
#   CPUS        taskset CPU list for the Loom tools               (default: all CPUs)
#   ROCM_PATH   ROCm install                                      (default /opt/rocm)
#   LLVM        dir with llvm-readelf   (default $ROCM_PATH/llvm/bin, else $ROCM_PATH/lib/llvm/bin)
#   IMPORT_TIMEOUT, COMPILE_TIMEOUT  seconds                      (default 1800, 3600)
#   PROBE       optional managed-memory probe binary to run first
#   GPU_WAIT    seconds to wait for other users' GPU processes    (default 900)
#   OUT         summary file                                      (default $TEST_DIR/compare.summary)
#   TMPDIR      kept if set; else $LOOM_REMOTE_BASE/tmp
#   LOOM_REMOTE_BASE  base dir                                    (default $HOME/loom-work)
#
# Output: $OUT (default $TEST_DIR/compare.summary; one block per step) and
# loom_root.report.json (the compile report) in TEST_DIR.
set -u
usage() { sed -n '2,/^set -u/p' "$0" | sed '$d'; exit "${1:-2}"; }
case "${1:-}" in -h|--help) usage 0;; '') ;; *) usage;; esac
for v in TEST_DIR KERNEL_SRC INPUTS TOOLS; do
    [ -n "${!v:-}" ] || { echo "set $v"; usage; }
done
ROOT=${ROOT:-mpl::burn_zones}
NCELL=${NCELL:-16}
RUNS=${RUNS:-}
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
GPU_WAIT=${GPU_WAIT:-900}
OUT=${OUT:-$(cd "$TEST_DIR" && pwd)/compare.summary}   # never next to this script: the skill dir may be read-only
export TMPDIR=${TMPDIR:-${LOOM_REMOTE_BASE:-$HOME/loom-work}/tmp}
mkdir -p "$TMPDIR"
say() { echo "$*" | tee -a "$OUT"; }
: > "$OUT"
say "# $(date) $(hostname) target=$TARGET ncell=$NCELL root=$ROOT"
cd "$TEST_DIR" || exit 1

# 1. Loom compile of the real root for NCELL^3 zones
ZONES=$((NCELL * NCELL * NCELL))
WG=$(( (ZONES + 127) / 128 ))
LROOT=${ROOT//::/.}
if [ -f "$TOOLS/tools.md5" ]; then (cd "$TOOLS" && md5sum -c --quiet tools.md5) || { say "TOOLS MD5 MISMATCH"; exit 1; }; fi
say "== loom compile: zones=$ZONES workgroup_count.x=$WG tools=$(md5sum "$TOOLS/loom-compile" | cut -c1-12)"
t0=$(date +%s)
if ! taskset -c "$CPUS" timeout "$IMPORT_TIMEOUT" "$TOOLS/loom-import-cxx" --data-model=lp64 --fp-model=precise \
      --root="$ROOT" --output=loom_root.loom "$KERNEL_SRC" > loom_import.log 2>&1; then
    say "IMPORT FAIL"; grep -m5 error loom_import.log | tee -a "$OUT"; exit 1
fi
say "import ok ($(( $(date +%s) - t0 ))s)"
t0=$(date +%s)
if ! taskset -c "$CPUS" timeout "$COMPILE_TIMEOUT" "$TOOLS/loom-compile" loom_root.loom --root="$LROOT" \
      --config="$LROOT.workgroup_count.x=$WG" --config="$LROOT.workgroup_count.y=1" --config="$LROOT.workgroup_count.z=1" \
      --format=amdgpu-hsaco --target=amdgpu:"$TARGET" --compile-report=summary \
      --compile-report-output=loom_root.report.json --output=loom_root.hsaco > loom_compile.log 2>&1; then
    say "COMPILE FAIL ($(( $(date +%s) - t0 ))s)"; grep -m5 error loom_compile.log | cut -c1-300 | tee -a "$OUT"; exit 1
fi
say "compile ok ($(( $(date +%s) - t0 ))s)"
"$LLVM/llvm-readelf" --notes loom_root.hsaco | awk '{for (j = 1; j < NF; j++) { if ($j == ".vgpr_count:") v = $(j+1); if ($j == ".agpr_count:") a = $(j+1); if ($j == ".private_segment_fixed_size:") s = $(j+1) }} END{printf "resources: vgpr=%s agpr=%s scratch=%sB\n", v, a, s}' | tee -a "$OUT"
say "symbols: $("$LLVM/llvm-readelf" --syms loom_root.hsaco | awk '$4 == "FUNC" || $4 == "OBJECT" {print $8}' | sort -u | tr '\n' ' ')"
say "kernel arguments (code-object metadata, in order):"
"$LLVM/llvm-readelf" --notes loom_root.hsaco | awk '/\.args:/{on=1} on && /\.name:|\.offset:|\.size:|\.value_kind:/{gsub(/^ +-? */,""); printf "%s ", $0} on && /\.value_kind:/{print ""} /\.group_segment_fixed_size|\.kernarg_segment_size/{on=0}' | tee -a "$OUT"

# 2. wait for other users' GPU processes (shared GPU; correctness only, but be polite)
# (no test_react process of ours is running yet, so every listed PID is someone else's)
waited=0 others="?"
SMI=$(command -v rocm-smi || ls "$ROCM_PATH/bin/rocm-smi" 2>/dev/null)
while [ -n "$SMI" ] && [ "$waited" -lt "$GPU_WAIT" ]; do
    others=$("$SMI" --showpids 2>/dev/null | awk '/^[0-9]+/{print $1}' | wc -l)
    [ "$others" -eq 0 ] && break
    sleep 30; waited=$((waited + 30))
done
[ -n "$SMI" ] || others="? (rocm-smi not found)"
say "== gpu: other processes after ${waited}s wait: $others"

# 3. optional managed-memory probe
if [ -n "${PROBE:-}" ] && [ -x "$PROBE" ]; then
    say "== managed-memory probe"; timeout 60 "$PROBE" 2>&1 | tee -a "$OUT"
fi

# 4. compare runs
EXE=$(ls -t ./*.hip*.ex 2>/dev/null | head -1)
[ -n "$EXE" ] || { say "no HIP executable in $TEST_DIR"; exit 1; }
IFS=';' read -r -a specs <<< "$RUNS"
[ ${#specs[@]} -gt 0 ] || specs=("")   # RUNS="" -> one run with the inputs file alone
n=0
for spec in "${specs[@]}"; do
    n=$((n + 1))
    log=compare_run$n.log
    say "== compare run $n: $EXE $INPUTS n_cell=$NCELL max_grid_size=$NCELL $spec"
    timeout 1800 "$EXE" "$INPUTS" n_cell="$NCELL" max_grid_size="$NCELL" integrator.use_jacobian_caching=0 \
        unit_test.loom_hsaco=loom_root.hsaco unit_test.loom_kernel="$LROOT" unit_test.loom_ncells="$ZONES" \
        unit_test.loom_compare=1 unit_test.loom_repeats=2 $spec > "$log" 2>&1
    rc=$?
    grep -E "^LOOM_COMPARE|^LOOM_TIMING|Abort|abort|error" "$log" | head -20 | tee -a "$OUT"
    say "exit code: $rc"
done
say "# done $(date)"
