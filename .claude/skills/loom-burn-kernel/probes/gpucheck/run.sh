#!/bin/bash
# GPU result check for probes: run each probe three ways and compare the 64
# outputs bit for bit. A probe that imports and compiles (probe.sh) can still
# be miscompiled; this checks that it computes what hipcc computes.
#   host   g++ -O0 -ffp-contract=off, kernel body run for threadIdx.x = 0..63
#          (fakehip/hip/hip_runtime.h stands in for HIP)
#   hipcc  ROCm clang, plain gfx90a code object; the probe's KERNEL macro is
#          swapped for __launch_bounds__(64) and `static void` helpers get
#          __device__ (Loom needs neither)
#   loom   loom-import-cxx + loom-compile, as probe.sh does
# Each code object runs as one 64-thread workgroup (gpu_driver.cpp) on inputs
# from common.H (negative, below and above 1). PASS = all three identical.
# Probes must have the kit's signature: k(double* out, const double* in), or
# k(double* out, const double* in, const int* sp).
#
# Usage (on the GPU host):
#   TOOLS=<dir with loom-import-cxx, loom-compile> bash run.sh <probe.cpp>...
# Writes out/ next to this script. ROCM_PATH defaults to /opt/rocm; the
# target is gfx90a (edit --offload-arch/--target for another GPU).
#
# Results at hrx-system main 920eb7f3ff (2026-11-07): d14 d22 f02 l01 o08 o09
# r08 r10 r12 r25 and the stricter variants l01b o08b o09b r08b r10b r12b all
# PASS (probes/README.md).
set -u
[ -n "${TOOLS:-}" ] || { sed -n '2,/^set -u/p' "$0" | sed '$d'; exit 2; }
G=$(cd "$(dirname "$0")" && pwd); O=$G/out; mkdir -p $O
L=${ROCM_PATH:-/opt/rocm}/lib/llvm/bin; [ -x $L/llvm-readelf ] || L=${ROCM_PATH:-/opt/rocm}/llvm/bin
hipcc -std=c++20 -O2 -o $O/gpu_driver $G/gpu_driver.cpp || exit 1
for src in "$@"; do n=$(basename ${src%.cpp}); src=$(realpath $src)
  nargs=2; grep -q "const int\* sp" $src && nargs=3
  if ! g++ -std=c++20 -O0 -ffp-contract=off -w -I $G/fakehip -include $src -DNARGS=$nargs -o $O/$n.host $G/host_driver.cpp 2> $O/$n.host.log; then
    printf '%-34s HOST-BUILD-FAIL %s\n' $n "$(grep -m1 error $O/$n.host.log)"; continue; fi
  $O/$n.host $O/$n.host.txt
  sed -e 's/^#define KERNEL .*/#define KERNEL __global__ __launch_bounds__(64)/' -e 's/^static void /static __device__ void /' $src > $O/$n.hip.cpp
  $L/clang++ -x hip -std=c++20 -O2 -ffp-contract=off --offload-arch=gfx90a --cuda-device-only --no-gpu-bundle-output -o $O/$n.hip.hsaco $O/$n.hip.cpp > $O/$n.hipcc.log 2>&1 || { printf '%-34s HIPCC-FAIL\n' $n; continue; }
  sym=$($L/llvm-readelf --syms $O/$n.hip.hsaco | awk '$4=="FUNC" && $8 ~ /5probe1k/ {print $8; exit}')
  $O/gpu_driver $O/$n.hip.hsaco "$sym" $nargs $O/$n.hip.txt || { printf '%-34s HIP-RUN-FAIL\n' $n; continue; }
  $TOOLS/loom-import-cxx --data-model=lp64 --fp-model=precise --root=probe::k --output=$O/$n.loom $src > $O/$n.import.log 2>&1 || { printf '%-34s LOOM-IMPORT-FAIL %s\n' $n "$(grep -m1 -o 'error.*' $O/$n.import.log | cut -c1-120)"; continue; }
  $TOOLS/loom-compile $O/$n.loom --root=probe.k --config=probe.k.workgroup_count.x=1 --config=probe.k.workgroup_count.y=1 \
     --config=probe.k.workgroup_count.z=1 --format=amdgpu-hsaco --target=amdgpu:gfx90a --output=$O/$n.loom.hsaco > $O/$n.compile.log 2>&1 \
     || { printf '%-34s LOOM-COMPILE-FAIL %s\n' $n "$(grep -m1 -o 'error.*' $O/$n.compile.log | cut -c1-120)"; continue; }
  $O/gpu_driver $O/$n.loom.hsaco probe.k $nargs $O/$n.loom.txt || { printf '%-34s LOOM-RUN-FAIL\n' $n; continue; }
  hh=$(cmp -s $O/$n.host.txt $O/$n.hip.txt && echo same || echo DIFF)
  hl=$(cmp -s $O/$n.hip.txt $O/$n.loom.txt && echo same || echo DIFF)
  ndiff=$(paste -d' ' <(awk '{print $2}' $O/$n.hip.txt) <(awk '{print $2}' $O/$n.loom.txt) | awk '$1!=$2' | wc -l)
  verdict=PASS; [ $hh = same ] && [ $hl = same ] || verdict=FAIL
  printf '%-34s %s  host~hipcc=%s hipcc~loom=%s (%d/64 lanes differ)  lane5=%s\n' $n $verdict $hh $hl $ndiff "$(awk '$1==5{print $3}' $O/$n.loom.txt)"
done
