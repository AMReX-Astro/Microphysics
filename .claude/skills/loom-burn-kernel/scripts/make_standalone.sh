#!/bin/bash
# Bundle a port's kernel tree into a self-contained standalone HIP reproducer
# (assets/standalone/README.md): no Microphysics, no AMReX, only hipcc.
#
#   make_standalone.sh [--tar] PORT_DIR OUT_DIR [CASE_FILE ...]
#
#   PORT_DIR   the port's kernel directory, unit_test/test_react/loom/<port>
#   OUT_DIR    created (must not exist yet); gets the kernel headers,
#              burn_kernel.cpp, main.hip, loom_case_format.H, Makefile,
#              README.md, shim/ (for MATH=recipes), cases/, MANIFEST
#   CASE_FILE  case files written by test_react (unit_test.loom_dump_case)
#              or by burn_standalone --out / --write-subset; copied to cases/
#   --tar      also write OUT_DIR.tar.gz
#
# Every #include in the bundle must resolve inside it (or be a HIP/C++
# standard header); the script stops otherwise. The validators
# (validate_*.H) are left out. The port must have the standalone mode: the
# MPL_HIP_STANDALONE branch in support.H and the MPL_GPU_KERNEL guards in
# burn_kernel.cpp (assets/kernel/support.H, assets/root/burn_kernel.cpp).
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
skill=$(cd "$here/.." && pwd)

tar=0
if [ "${1:-}" = "--tar" ]; then tar=1; shift; fi
if [ $# -lt 2 ]; then
  sed -n '2,19p' "$0" | sed 's/^# \{0,1\}//' >&2
  exit 2
fi
port=$(cd "$1" && pwd)
out=$2
shift 2

die() { echo "make_standalone.sh: $*" >&2; exit 1; }

[ -f "$port/burn_kernel.cpp" ] || die "$port has no burn_kernel.cpp (not a port kernel directory?)"
[ -f "$port/support.H" ] || die "$port has no support.H"
grep -q 'MPL_HIP_STANDALONE' "$port/support.H" ||
  die "$port/support.H has no MPL_HIP_STANDALONE branch: the port predates the standalone mode. Bring over the MPL_HIP_STANDALONE / MPL_HIP_MATH_RECIPES / MPL_GPU_KERNEL parts of $skill/assets/kernel/support.H (diff the two files)."
grep -q 'MPL_GPU_KERNEL' "$port/burn_kernel.cpp" ||
  die "$port/burn_kernel.cpp still selects the GPU zone index with __loom__: use MPL_GPU_KERNEL as in $skill/assets/root/burn_kernel.cpp."
[ ! -e "$out" ] || die "$out already exists"

mkdir -p "$out/shim/loomcxx" "$out/cases"
out=$(cd "$out" && pwd)

# kernel tree: every header except the validators, and the root
for f in "$port"/*.H; do
  case $(basename "$f") in validate_*) continue ;; esac
  cp "$f" "$out/"
done
cp "$port/burn_kernel.cpp" "$out/"
# driver, format, build, docs
cp "$skill/assets/standalone/main.hip" "$skill/assets/standalone/loom_case_format.H" \
   "$skill/assets/standalone/Makefile" "$skill/assets/standalone/README.md" "$out/"
cp "$skill/assets/f64-tests/ocml_gpu/device_shim/loomcxx/scalar.h" "$out/shim/loomcxx/"
for c in "$@"; do
  [ -f "$c" ] || die "case file $c not found"
  cp "$c" "$out/cases/"
done

# self-containment: every "..." include exists in the bundle, every <...>
# include is HIP, the C++ standard library, or the recipes' shim
allowed='^(hip/hip_runtime\.h|loomcxx/scalar\.h|algorithm|array|cinttypes|cmath|cstdarg|cstddef|cstdint|cstdio|cstdlib|cstring|limits|string|type_traits|utility|vector)$'
bad=0
while IFS= read -r line; do
  file=${line%%:*}
  inc=$(printf '%s\n' "$line" | sed -E 's/^[^:]*:[[:space:]]*#[[:space:]]*include[[:space:]]*//')
  case $inc in
    \"*\")
      name=${inc#\"}; name=${name%\"}
      [ -f "$out/$name" ] || { echo "  $(basename "$file"): $inc not in the bundle" >&2; bad=1; } ;;
    \<*\>)
      name=${inc#<}; name=${name%>}
      printf '%s\n' "$name" | grep -Eq "$allowed" || { echo "  $(basename "$file"): $inc is not a HIP/standard header" >&2; bad=1; } ;;
    *) echo "  $(basename "$file"): computed include $inc" >&2; bad=1 ;;
  esac
done < <(grep -H -E '^[[:space:]]*#[[:space:]]*include' "$out"/*.H "$out"/*.cpp "$out"/main.hip)
[ $bad -eq 0 ] || die "the bundle is not self-contained (see above); $out left for inspection"

{
  echo "standalone burn reproducer bundle"
  echo "created: $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
  echo "port: $port"
  echo "port git: $(git -C "$port" rev-parse HEAD 2>/dev/null || echo unknown)$(git -C "$port" status --porcelain -- . 2>/dev/null | grep -q . && echo ' (port files modified or untracked)')"
  echo "skill: $skill"
  echo "skill git: $(git -C "$skill" rev-parse HEAD 2>/dev/null || echo unknown)$(git -C "$skill" status --porcelain -- . 2>/dev/null | grep -q . && echo ' (skill files modified or untracked)')"
  echo
  echo "md5:"
  (cd "$out" && find . -type f ! -name MANIFEST | LC_ALL=C sort | while read -r f; do
     if command -v md5sum >/dev/null; then md5sum "$f"; else md5 -r "$f"; fi
   done)
} > "$out/MANIFEST"

echo "bundle: $out ($(ls "$out" | wc -l | tr -d ' ') entries, $(ls "$out/cases" | wc -l | tr -d ' ') case files)"
if [ $tar -eq 1 ]; then
  tar -C "$(dirname "$out")" -czf "$out.tar.gz" "$(basename "$out")"
  echo "tarball: $out.tar.gz"
fi
echo "build: cd $out && make [ARCH=gfxNNN] [MATH=recipes]"
echo "run:   ./burn_standalone_ocml cases/<case>"
