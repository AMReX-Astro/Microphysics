#!/bin/bash
# Build and run the GPU bit-for-bit check of the transcribed OCML f64 recipes
# against hipcc's OCML (see README.md).
#
# Runs ON the GPU host (needs hipcc, from PATH or $ROCM_PATH/bin, default
# /opt/rocm/bin, and a free GPU). Builds against
# the skill's own recipe headers in ../../kernel/.
#
# Usage: ./build_and_run.sh [samples_per_sweep]     (default 1048576)
#        ./build_and_run.sh -h                       (this text)
# Environment:
#   BUILD_DIR               output directory (default ${TMPDIR:-/tmp}/loom_ocml_gpu)
#   OFFLOAD_ARCH            default: first gfx* target from amdgpu-arch or
#                           rocminfo, else gfx90a
#   OCMLGPU_APPROX_VARIANT  lowering of the f32 approximate log2/exp2 used by
#                           cbrt (0 = __builtin_amdgcn_logf/exp2f, default;
#                           1 = __ocml_native_*; 2 = __builtin_log2f/exp2f)
#   RUN_TIMEOUT             seconds (default 300; GPU time itself is well
#                           under a second, most of the run is host-side
#                           input generation and comparison)
set -euo pipefail
case "${1:-}" in -h|--help) sed -n '2,/^set -euo/p' "$0" | sed '$d'; exit 0;; esac
here=$(cd "$(dirname "$0")" && pwd)
recipes=$(cd "$here/../../kernel" && pwd)
out=${BUILD_DIR:-${TMPDIR:-/tmp}/loom_ocml_gpu}
rocm=${ROCM_PATH:-/opt/rocm}
if [ -z "${OFFLOAD_ARCH:-}" ]; then
    OFFLOAD_ARCH=$( { "$rocm/llvm/bin/amdgpu-arch" 2>/dev/null || amdgpu-arch 2>/dev/null ||
                      "$rocm/bin/rocminfo" 2>/dev/null || rocminfo 2>/dev/null; } |
                    grep -m1 -o 'gfx[0-9a-f]*' || true)
fi
arch=${OFFLOAD_ARCH:-gfx90a}
variant=${OCMLGPU_APPROX_VARIANT:-0}
n=${1:-1048576}
mkdir -p "$out"

if ! command -v hipcc >/dev/null && [ -x "$rocm/bin/hipcc" ]; then PATH=$rocm/bin:$PATH; fi
command -v hipcc >/dev/null || { echo "hipcc not found in PATH or $rocm/bin (set ROCM_PATH)" >&2; exit 2; }
echo "offload arch: $arch"
echo "hipcc: $(command -v hipcc)"
hipcc --version 2>&1 | head -4 || true

# The stock Microphysics HIP reference flags; OCML itself is prebuilt bitcode
# and is not affected by -ffp-contract.
flags=(-std=c++20 -O3 -ffp-contract=off --offload-arch="$arch"
       -DOCMLGPU_APPROX_VARIANT="$variant"
       -I"$here/device_shim" -I"$recipes")

hipcc "${flags[@]}" "$here/ocml_gpu_check.hip" -o "$out/ocml_gpu_check"

# Device assembly, for inspecting what each side lowered to (optional).
asm="$out/ocml_gpu_check.$arch.s"
if hipcc "${flags[@]}" --offload-device-only -S "$here/ocml_gpu_check.hip" -o "$asm" 2>/dev/null ||
   hipcc "${flags[@]}" --cuda-device-only -S "$here/ocml_gpu_check.hip" -o "$asm" 2>/dev/null; then
    echo "device asm: $asm"
    for ins in v_log_f32 v_exp_f32 v_rcp_f64 v_div_scale_f64 v_trunc_f64 v_rndne_f64 \
               v_ldexp_f64 'swappc'; do
        printf '  %-16s %s\n' "$ins" "$(grep -c "$ins" "$asm" || true)"
    done
else
    echo "(device asm dump failed; ignored)"
fi

set +e
timeout "${RUN_TIMEOUT:-300}" "$out/ocml_gpu_check" "$n" | tee "$out/ocml_gpu_check.out"
status=${PIPESTATUS[0]}
echo "output saved to $out/ocml_gpu_check.out (exit status $status)"
exit "$status"
