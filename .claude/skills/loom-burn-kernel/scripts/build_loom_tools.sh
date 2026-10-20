#!/bin/bash
# Build a frozen loom-import-cxx + loom-compile pair from an hrx-system source
# tree, optionally with local patches (see ../patches/README.md).
#
# Runs on the GPU host (locally, or there over one ssh connection), detached;
# keep the backgrounded command on its own line, stdin from /dev/null:
#   SRC=... OUT=... PATCHES="a.patch b.patch" \
#   nohup setsid bash build_loom_tools.sh > build.log 2>&1 < /dev/null &
# The full build takes tens of minutes; check build.log or the BUILDING marker.
#
# Steps: copy SRC to $OUT/src (the original tree is never touched; .git and
# ./build, ./build-* are excluded), apply each patch with `patch -p1`, configure with
# CMake/Ninja using ROCm's clang and a C++23 libstdc++ (GCC >= 14), build the two tools,
# copy them to $OUT/tools and write $OUT/tools/tools.md5 (checked by
# loom_smoke.sh, loom_gpu_compare.sh, probe.sh) and $OUT/tools/SOURCE (source
# commit, dirty state, patches and their md5).
#
# Usage: build_loom_tools.sh            (all configuration via environment)
#        build_loom_tools.sh -h         (this text)
#
# Fresh machine (Ubuntu 24.04 + ROCm): apt install g++-14 libstdc++-14-dev
# python3.12-venv patch git, then
#   git clone https://github.com/ROCm/hrx-system.git && cd hrx-system
#   git checkout fada6788a1      # pinned compiler (SKILL.md, Phase 1)
#   python3 dev.py cmake setup   # .venv with the pinned cmake/ninja (used below)
#
# Environment:
#   LOOM_REMOTE_BASE  base dir for the defaults below      (default $LOOM_WORK, else $HOME/loom-work)
#   SRC            hrx-system source tree                 (default $LOOM_REMOTE_BASE/hrx-system)
#   OUT            work dir; src/, build/, tools/ inside   (default $LOOM_REMOTE_BASE/loom-tools/<UTC timestamp>)
#   PATCHES        space-separated patch files, applied in order with patch -p1 (default none)
#   ROCM           ROCm install with llvm/bin/clang        (default $ROCM_PATH, else /opt/rocm)
#   PYVENV         bin dir with cmake, ninja, python3      (default $SRC/.venv/bin from `dev.py cmake
#                                                          setup` if present, else the tools on PATH)
#   GCC_TOOLCHAIN  GCC >= 14 prefix (lib/gcc/<triple>/<ver>/), e.g. a conda gcc14; its lib/ is
#                  added as rpath                          (default: none, clang picks the newest
#                                                          system GCC, e.g. Ubuntu's g++-14)
#   GCC_INSTALL_DIR  --gcc-install-dir for clang           (default: glob $GCC_TOOLCHAIN/lib/gcc/*/*/, newest)
#   TARGETS        IREE_HAL_AMDGPU_TARGETS                 (default: the local GPU from amdgpu-arch /
#                                                          rocminfo, else gfx90a)
#   JOBS           build parallelism                       (default nproc)
#   CPUS           taskset CPU list for the build          (default all CPUs)
#   TMPDIR         temporary dir                           (default $OUT/tmp)
#
# moth: LOOM_REMOTE_BASE=/mnt/data/bwibking ROCM=/opt/rocm/core-10.1
#       PYVENV=$LOOM_REMOTE_BASE/venv-py312/bin GCC_TOOLCHAIN=$LOOM_REMOTE_BASE/gcc14
set -u
case "${1:-}" in -h|--help) sed -n '2,/^set -u/p' "$0" | sed '$d'; exit 0;; '') ;; *) sed -n '2,/^set -u/p' "$0" | sed '$d'; exit 2;; esac

B=${LOOM_REMOTE_BASE:-${LOOM_WORK:-$HOME/loom-work}}
SRC=${SRC:-$B/hrx-system}
OUT=${OUT:-$B/loom-tools/$(date -u +%Y%m%dT%H%M%SZ)}
PATCHES=${PATCHES:-}
ROCM=${ROCM:-${ROCM_PATH:-/opt/rocm}}
GCC14=${GCC_TOOLCHAIN:-}
detect_gfx() {
    { "$ROCM/llvm/bin/amdgpu-arch" || "$ROCM/bin/rocminfo"; } 2>/dev/null | grep -m1 -o 'gfx[0-9a-f]*'
}
TARGETS=${TARGETS:-$(detect_gfx)}
[ -n "$TARGETS" ] || { echo "WARNING: no TARGETS given and no GPU detected; building for gfx90a" >&2; TARGETS=gfx90a; }
JOBS=${JOBS:-$(nproc)}
CPUS=${CPUS:-0-$(( $(nproc) - 1 ))}

T0=$(date +%s)
log() { echo "[$(date '+%F %T') +$(( $(date +%s) - T0 ))s] $*"; }
fail() { log "FAILED: $*"; exit 1; }

[ -d "$SRC" ] || fail "SRC=$SRC is not a directory (git clone https://github.com/ROCm/hrx-system.git; see -h)"
if [ -n "${PYVENV:-}" ]; then V=$PYVENV
elif [ -x "$SRC/.venv/bin/cmake" ]; then V=$SRC/.venv/bin
else V=$(dirname "$(command -v cmake || echo /usr/bin/cmake)"); fi
PY=$(command -v "$V/python3.12" "$V/python3" python3.12 python3 2>/dev/null | head -1)
NINJA=$(command -v "$V/ninja" ninja 2>/dev/null | head -1)
for x in "$ROCM/llvm/bin/clang" "$ROCM/llvm/bin/clang++" "$V/cmake" "$NINJA" "$PY"; do
    [ -n "$x" ] && [ -x "$x" ] || fail "missing executable ${x:-ninja/python3} (check ROCM / PYVENV; run dev.py cmake setup in SRC)"
done
if [ -z "${GCC_INSTALL_DIR:-}" ] && [ -n "$GCC14" ]; then
    GCC_INSTALL_DIR=$(ls -d "$GCC14"/lib/gcc/*/*/ 2>/dev/null | sort -V | tail -1)
    GCC_INSTALL_DIR=${GCC_INSTALL_DIR%/}
    [ -n "$GCC_INSTALL_DIR" ] && [ -d "$GCC_INSTALL_DIR" ] || fail "no gcc install dir under $GCC14/lib/gcc/*/*/ (set GCC_TOOLCHAIN or GCC_INSTALL_DIR)"
fi
GCC_FLAG=${GCC_INSTALL_DIR:+--gcc-install-dir=$GCC_INSTALL_DIR}
RPATH=${GCC14:+-Wl,-rpath,$GCC14/lib}
# loom-import-cxx needs C++23 std::expected; check the toolchain before a long build
echo '#include <expected>
int main() { std::expected<int, int> e = 1; return *e - 1; }' > "${TMPDIR:-/tmp}/loom_expected_$$.cpp"
"$ROCM/llvm/bin/clang++" -std=c++23 $GCC_FLAG -fsyntax-only "${TMPDIR:-/tmp}/loom_expected_$$.cpp" 2>/dev/null \
    || { rm -f "${TMPDIR:-/tmp}/loom_expected_$$.cpp"; fail "ROCm clang++ finds no C++23 libstdc++ with std::expected (apt install g++-14 libstdc++-14-dev, or set GCC_TOOLCHAIN)"; }
rm -f "${TMPDIR:-/tmp}/loom_expected_$$.cpp"
abs_patches=""
for p in $PATCHES; do
    [ -f "$p" ] || fail "patch not found: $p"
    abs_patches="$abs_patches $(cd "$(dirname "$p")" && pwd)/$(basename "$p")"
done
PATCHES=${abs_patches# }

mkdir -p "$OUT" || fail "cannot create OUT=$OUT"
OUT=$(cd "$OUT" && pwd)
if [ -e "$OUT/BUILDING" ] && kill -0 "$(cat "$OUT/BUILDING" 2>/dev/null)" 2>/dev/null; then
    echo "BUSY: $OUT has a running build (PID $(cat "$OUT/BUILDING"))"; exit 3
fi
echo $$ > "$OUT/BUILDING"
trap 'rm -f "$OUT/BUILDING"' EXIT
export TMPDIR=${TMPDIR:-$OUT/tmp}
mkdir -p "$TMPDIR"

BUILD=$OUT/build
TOOLS=$OUT/tools
log "host=$(hostname) SRC=$SRC OUT=$OUT targets=$TARGETS jobs=$JOBS cpus=$CPUS"
log "ROCM=$ROCM PYVENV=$V GCC_INSTALL_DIR=${GCC_INSTALL_DIR:-<system>}"

# Source provenance (before copying, from the original tree).
commit=unknown dirty=""
if git -C "$SRC" rev-parse HEAD > /dev/null 2>&1; then
    commit=$(git -C "$SRC" rev-parse HEAD)
    [ -n "$(git -C "$SRC" status --porcelain --untracked-files=no 2>/dev/null)" ] && dirty=" (tracked files modified)"
elif [ -f "$SRC/COMMIT" ]; then
    commit=$(cat "$SRC/COMMIT")
fi
log "source commit: $commit$dirty"

# 1. fresh copy of the source tree
log "copying source -> $OUT/src"
rm -rf "$OUT/src"
mkdir -p "$OUT/src"
(cd "$SRC" && tar --exclude=./.git --exclude=./.venv --exclude=./build --exclude='./build-*' -cf - .) | (cd "$OUT/src" && tar -xf -) \
    || fail "copying $SRC"
echo "$commit" > "$OUT/src/COMMIT"
log "copy done"

# 2. patches
for p in $PATCHES; do
    log "applying $p"
    (cd "$OUT/src" && patch -p1 --forward --batch < "$p") > "$OUT/patch.$(basename "$p").log" 2>&1 \
        || { cat "$OUT/patch.$(basename "$p").log"; fail "patch $p does not apply to $commit"; }
done

# 3. configure
log "configuring"
"$V/cmake" -G Ninja -S "$OUT/src" -B "$BUILD" \
    -DCMAKE_MAKE_PROGRAM="$NINJA" -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER="$ROCM/llvm/bin/clang" -DCMAKE_CXX_COMPILER="$ROCM/llvm/bin/clang++" \
    "-DCMAKE_C_FLAGS=$GCC_FLAG" "-DCMAKE_CXX_FLAGS=$GCC_FLAG" \
    "-DCMAKE_EXE_LINKER_FLAGS=$RPATH" "-DCMAKE_SHARED_LINKER_FLAGS=$RPATH" \
    "-DCMAKE_MODULE_LINKER_FLAGS=$RPATH" \
    -DPython3_EXECUTABLE="$PY" -DIREE_ROCM_PATH="$ROCM" \
    -DIREE_HAL_DRIVER_DEFAULTS=OFF -DIREE_HAL_DRIVER_AMDGPU=ON -DIREE_HAL_DRIVER_TASK=ON \
    -DIREE_HAL_AMDGPU_TARGETS="$TARGETS" -DLOOM_IMPORT_CXX=ON \
    > "$OUT/configure.log" 2>&1 \
    || { tail -30 "$OUT/configure.log"; fail "configure (full log: $OUT/configure.log)"; }
log "configure done"

# 4. build
log "building loom-import-cxx loom-compile"
PATH=$V:$(dirname "$NINJA"):$PATH nice -n 10 taskset -c "$CPUS" "$V/cmake" --build "$BUILD" \
    --target loom-import-cxx loom-compile -j "$JOBS" > "$OUT/build.ninja.log" 2>&1 \
    || { grep -m20 -E 'error|FAILED' "$OUT/build.ninja.log"; fail "build (full log: $OUT/build.ninja.log)"; }
log "build done"

# 5. install the frozen pair
mkdir -p "$TOOLS"
cp "$BUILD/loom/src/loom/tools/loom-import-cxx/loom-import-cxx" \
   "$BUILD/loom/src/loom/tools/loom-compile/loom-compile" "$TOOLS/" || fail "binaries not found under $BUILD/loom/src/loom/tools"
(cd "$TOOLS" && md5sum loom-import-cxx loom-compile > tools.md5) || fail "md5sum"
{
    echo "source: $SRC"
    echo "commit: $commit$dirty"
    echo "built: $(date) on $(hostname)"
    echo "targets: $TARGETS"
    for p in $PATCHES; do echo "patch: $(md5sum "$p")"; done
} > "$TOOLS/SOURCE"
log "tools: $TOOLS"
sed 's/^/  /' "$TOOLS/tools.md5"
log "DONE"
