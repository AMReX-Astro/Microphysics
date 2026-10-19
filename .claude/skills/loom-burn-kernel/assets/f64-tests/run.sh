#!/bin/sh
# Build and run the host validation of f64_math_extra.H (see README.md).
# Usage: ./run.sh [samples_per_sweep]   (default 1000000)
# Needs clang++ (CXX overrides; falls back to ROCm's) and MPFR/GMP
# (macOS: brew install mpfr; Ubuntu: apt install libmpfr-dev libgmp-dev clang).
set -e
here=$(cd "$(dirname "$0")" && pwd)
recipes="$here/../kernel"
out=${BUILD_DIR:-${TMPDIR:-/tmp}/loom_ocml_test}
mkdir -p "$out"
mpfr=${MPFR_PREFIX:-$(brew --prefix mpfr 2>/dev/null || echo /usr)}
gmp=${GMP_PREFIX:-$(brew --prefix gmp 2>/dev/null || echo /usr)}
cxx=${CXX:-clang++}
if ! command -v "$cxx" >/dev/null 2>&1 && [ -x "${ROCM_PATH:-/opt/rocm}/llvm/bin/clang++" ]; then
  cxx=${ROCM_PATH:-/opt/rocm}/llvm/bin/clang++
fi
# system prefixes (/usr) are already on the default search paths
inc="" lib=""
for p in "$mpfr" "$gmp"; do
  [ "$p" = /usr ] && continue
  inc="$inc -I$p/include"; lib="$lib -L$p/lib"
done
for opt in O2 O0; do
  # shellcheck disable=SC2086  # $inc/$lib word-split on purpose
  "$cxx" -std=c++20 -$opt -ffp-contract=off -Wall -Wextra \
    -I"$here/host_shim" -I"$recipes" $inc \
    "$here/test_ocml.cpp" $lib -lmpfr -lgmp \
    -o "$out/test_ocml_$opt"
done
"$out/test_ocml_O2" "${1:-1000000}"
# The -O0 build must produce the same bits (output hash) as -O2.
"$out/test_ocml_O0" 20000 | grep 'output hash' > "$out/h0"
"$out/test_ocml_O2" 20000 | grep 'output hash' > "$out/h2"
cmp "$out/h0" "$out/h2" && echo "O0/O2 output hashes agree: $(cat "$out/h2")"
