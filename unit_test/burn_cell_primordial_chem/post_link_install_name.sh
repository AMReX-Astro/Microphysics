#!/bin/sh
#
# Wrapper around the AMReX linker command that patches stale Zerobrew GCC
# install_names after linking on macOS.

set -eu

real_linker="$1"
shift

output=""
prev=""
for arg in "$@"; do
    if [ "$prev" = "-o" ]; then
        output="$arg"
        break
    fi
    prev="$arg"
done

"$real_linker" "$@"

if [ "$(uname)" != "Darwin" ]; then
    exit 0
fi

if [ -z "$output" ] || [ ! -f "$output" ]; then
    exit 0
fi

if ! command -v otool >/dev/null 2>&1; then
    exit 0
fi

if ! command -v install_name_tool >/dev/null 2>&1; then
    exit 0
fi

target_root="/opt/zerobrew/prefix/lib/gcc/current"

for lib in libgfortran.5.dylib libquadmath.0.dylib libstdc++.6.dylib; do
    new_path="${target_root}/${lib}"
    if [ ! -e "$new_path" ]; then
        continue
    fi

    old_path=$(otool -L "$output" | awk -v lib="$lib" 'index($1, lib) {print $1; exit}')
    if [ -n "${old_path}" ] && [ "${old_path}" != "${new_path}" ]; then
        install_name_tool -change "$old_path" "$new_path" "$output"
    fi
done
