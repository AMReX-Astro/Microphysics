#!/bin/bash
# Poor man's profiler for a long-running loom-compile: sample its call stack
# with eu-stack every INTERVAL seconds until it exits, to see which pass
# (allocation, repair, scheduling, ...) a slow compile is stuck in.
#
# Runs on the machine running loom-compile (the local GPU box, or the remote
# host via one ssh connection; background only the job, on its own line):
#   cd <samples dir> || exit 1
#   nohup setsid bash loom_stack_sampler.sh [PID] > sampler.log 2>&1 < /dev/null &
#
# Needs eu-stack (Ubuntu: apt install elfutils) and permission to ptrace the
# compile. Ubuntu's default kernel.yama.ptrace_scope=1 forbids attaching to a
# process that isn't a child of the sampler, so every sample fails unless
# ptrace_scope is 0 (sudo sysctl -w kernel.yama.ptrace_scope=0) or the
# sampler runs as root; the script warns at start and logs eu-stack's error.
#
# Usage: loom_stack_sampler.sh [PID]
#        loom_stack_sampler.sh -h          (this text)
#   PID: the loom-compile process. Without it, the one loom-compile owned by
#   this user is used (pgrep -u $(id -u) -x loom-compile, optionally only
#   those whose working directory is CWD); more than one match is an error.
#   pgrep -x matches the exact process name, so it never matches this script.
#   Never kill compiles by pattern (pkill -f): other users share the host.
#
# Output (in OUT): s.NNN per sample (date, ps etime/time/rss, then the
# innermost DEPTH frames of every thread, one function per line), and
# samples.log with one line per sample: NNN time etime cputime rss and the
# innermost SUMMARY_FRAMES frames of the first thread joined by " < ".
#
# Environment:
#   INTERVAL        seconds between samples                  (default 60)
#   DEPTH           frames kept per sample (head -N)         (default 40)
#   SUMMARY_FRAMES  frames in the samples.log line           (default 4)
#   OUT             output dir                               (default .)
#   CWD             with no PID: only a loom-compile whose cwd is this dir
#   MAX_SAMPLES     stop after this many samples             (default unlimited)
set -u
usage() { sed -n '2,/^set -u/p' "$0" | sed '$d'; exit "${1:-2}"; }
case "${1:-}" in -h|--help) usage 0;; esac
[ $# -le 1 ] || usage

INTERVAL=${INTERVAL:-60}
DEPTH=${DEPTH:-40}
SUMMARY_FRAMES=${SUMMARY_FRAMES:-4}
OUT=${OUT:-.}
MAX_SAMPLES=${MAX_SAMPLES:-0}
command -v eu-stack > /dev/null || { echo "eu-stack (elfutils) not found in PATH (Ubuntu: apt install elfutils)"; exit 2; }
ptrace_scope=$(cat /proc/sys/kernel/yama/ptrace_scope 2>/dev/null || echo 0)
if [ "$ptrace_scope" -gt 0 ] && [ "$(id -u)" -ne 0 ]; then
    echo "WARNING: kernel.yama.ptrace_scope=$ptrace_scope: eu-stack -p cannot attach to a non-child process"
    echo "  as a normal user; samples will fail. Fix: sudo sysctl -w kernel.yama.ptrace_scope=0, or run with sudo."
fi

if [ $# -eq 1 ]; then
    P=$1
    case "$P" in *[!0-9]*|'') echo "PID must be a number"; usage;; esac
else
    cands=()
    for p in $(pgrep -u "$(id -u)" -x loom-compile); do
        if [ -n "${CWD:-}" ]; then
            [ "$(readlink "/proc/$p/cwd" 2>/dev/null)" = "$(cd "$CWD" 2>/dev/null && pwd -P)" ] || continue
        fi
        cands+=("$p")
    done
    if [ ${#cands[@]} -ne 1 ]; then
        echo "found ${#cands[@]} loom-compile processes of $(id -un)${CWD:+ in $CWD}: ${cands[*]:-none}; pass a PID"
        exit 2
    fi
    P=${cands[0]}
fi
kill -0 "$P" 2>/dev/null || { echo "process $P not running (or not ours)"; exit 2; }

mkdir -p "$OUT" || exit 1
cd "$OUT" || exit 1
echo "# $(date) sampling PID $P ($(ps -o comm= -p "$P")) every ${INTERVAL}s; cwd $(readlink "/proc/$P/cwd" 2>/dev/null)" | tee -a samples.log
n=0
while kill -0 "$P" 2>/dev/null; do
    n=$((n + 1))
    f=$(printf 's.%03d' "$n")
    ps_line=$(ps -o etime=,time=,rss= -p "$P" 2>/dev/null | awk '{print $1, $2, $3 "kB"}')
    {
        echo "# $(date '+%F %T') etime/time/rss: $ps_line"
        eu-stack -p "$P" 2> "$f.err" | awk '/^#/{print $NF}' | head -n "$DEPTH"
    } > "$f"
    frames=$(sed -n "2,$((SUMMARY_FRAMES + 1))p" "$f" | paste -sd'<' - | sed 's/</ < /g')
    if [ -z "$frames" ]; then
        frames="(no stack: $(head -1 "$f.err" | cut -c1-160))"
        [ -s "$f.err" ] || frames="(no stack, no eu-stack error)"
    fi
    [ -s "$f.err" ] || rm -f "$f.err"
    echo "$(printf '%03d' "$n") $(date '+%T') $ps_line $frames" >> samples.log
    [ "$MAX_SAMPLES" -gt 0 ] && [ "$n" -ge "$MAX_SAMPLES" ] && break
    sleep "$INTERVAL"
done
echo "# $(date) stopped after $n samples ($(kill -0 "$P" 2>/dev/null && echo 'MAX_SAMPLES reached' || echo 'process exited'))" | tee -a samples.log
