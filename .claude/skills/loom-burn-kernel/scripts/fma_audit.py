#!/usr/bin/env python3
"""Screen a Loom code object for FP contraction by accounting for its
v_fma_f64 count (SKILL.md Phase 6).

This is a heuristic, not an exact per-instruction attribution: the sqrt
remainder is accepted up to --max-per-sqrt per sqrt, and a negative
remainder is accepted as constant folding, so a few contracted operations
could hide in that slack. Its "OK" means "no sign of contraction"; the
`m_fp_contract` probe (probes/) is the direct evidence.

A raw v_fma_f64 count says nothing about FP contraction: AMDGPU lowers each
f64 division to 5 v_fma_f64 (next to v_div_scale/v_div_fmas/v_div_fixup) and
each f64 sqrt (v_rsq_f64 + Newton steps) to a few more, on top of the math
recipes' explicit fma calls. This script counts

  explicit  scalar.fmaf ops a fully inlined kernel contains: in the .loom IR,
            the kernel body plus every func.def reachable from it, summed over
            call paths (static count; loops are not unrolled)
  division  5 x v_div_fixup_f64 in the disassembly
  sqrt      the rest, compared against v_rsq_f64 / v_sqrt_f64

and reports the unexplained remainder per sqrt. Measured on triple_alpha
(hrx-system 4031286f0b, gfx90a): 32,866 = 19,401 explicit + 12,475 division
+ 990 for 303 sqrts (3.3 per sqrt), i.e. no contraction.

  fma_audit.py --loom loom_root.loom --hsaco loom_root.hsaco [--llvm DIR]

Exit status 0 if everything is accounted for (a small negative remainder is
constant-folded recipe calls: contraction only ever adds FMAs), 1 if the
remainder is above --max-per-sqrt per sqrt (default 8; then look for
contraction), 2 on usage errors.
"""
import argparse
import functools
import glob
import os
import re
import shutil
import subprocess
import sys

FMAS_PER_DIV = 5


def find_objdump(llvm):
    cands = []
    if llvm:
        cands.append(os.path.join(llvm, 'llvm-objdump'))
    rocm = os.environ.get('ROCM_PATH', '/opt/rocm')
    cands += [os.path.join(rocm, 'llvm/bin/llvm-objdump'), os.path.join(rocm, 'lib/llvm/bin/llvm-objdump')]
    cands += sorted(glob.glob('/opt/rocm*/llvm/bin/llvm-objdump')) + sorted(glob.glob('/opt/rocm*/lib/llvm/bin/llvm-objdump'))
    w = shutil.which('llvm-objdump')
    if w:
        cands.append(w)
    for c in cands:
        if os.access(c, os.X_OK):
            return c
    sys.exit('no llvm-objdump found (pass --llvm DIR or set ROCM_PATH)')


def explicit_fmas(loom_path):
    funcs = {'__kernel__': {'fma': 0, 'calls': []}}
    cur = '__kernel__'   # the kernel.def body precedes the first func.def
    with open(loom_path) as f:
        for line in f:
            m = re.match(r'\s*func\.def \w+ @(\S+?)\(', line)
            if m:
                cur = m.group(1)
                funcs[cur] = {'fma': 0, 'calls': []}
                continue
            if 'scalar.fmaf' in line:
                funcs[cur]['fma'] += 1
            funcs[cur]['calls'] += re.findall(r'func\.call(?: [a-z]+)* @([A-Za-z0-9_$.]+)\(', line)

    sys.setrecursionlimit(100000)

    @functools.lru_cache(None)
    def n(fn):
        d = funcs.get(fn)
        return 0 if d is None else d['fma'] + sum(n(c) for c in d['calls'])

    return n('__kernel__'), len(funcs) - 1


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--loom', required=True, help='the .loom IR the code object was compiled from')
    ap.add_argument('--hsaco', required=True, help='the Loom code object')
    ap.add_argument('--llvm', help='dir with llvm-objdump (default: ROCm, then PATH)')
    ap.add_argument('--max-per-sqrt', type=float, default=8.0,
                    help='largest plausible remainder per sqrt before flagging contraction (default 8)')
    a = ap.parse_args()
    for f in (a.loom, a.hsaco):
        if not os.path.isfile(f):
            ap.error(f'not a file: {f}')

    explicit, nfuncs = explicit_fmas(a.loom)
    dis = subprocess.run([find_objdump(a.llvm), '-d', a.hsaco], capture_output=True, text=True)
    if dis.returncode != 0:
        sys.exit(f'llvm-objdump failed: {dis.stderr.strip()[:300]}')
    # mnemonics may carry an encoding suffix (v_rsq_f64_e32)
    count = {k: len(re.findall(r'\b' + k + r'(?:_e32|_e64)?\b', dis.stdout))
             for k in ('v_fma_f64', 'v_div_fixup_f64', 'v_rsq_f64', 'v_sqrt_f64', 'v_mul_f64', 'v_add_f64')}
    division = FMAS_PER_DIV * count['v_div_fixup_f64']
    nsqrt = count['v_rsq_f64'] + count['v_sqrt_f64']
    rest = count['v_fma_f64'] - explicit - division

    print(f'v_fma_f64 total            {count["v_fma_f64"]}')
    print(f'  explicit fma (IR paths)  {explicit}   ({nfuncs} functions in the IR)')
    print(f'  f64 division (5 x {count["v_div_fixup_f64"]})  {division}')
    print(f'  remainder                {rest}   vs {nsqrt} sqrt ({count["v_rsq_f64"]} v_rsq_f64, {count["v_sqrt_f64"]} v_sqrt_f64)')
    print(f'(v_mul_f64 {count["v_mul_f64"]}, v_add_f64 {count["v_add_f64"]})')
    per = rest / nsqrt if nsqrt else float('inf') if rest > 0 else 0.0
    if rest < 0:
        # contraction only adds FMAs; fewer means recipe calls with constant
        # arguments were folded away (seen in smoke_rates: exp(const))
        print(f'OK: {-rest} fewer FMAs than explicit + division (constant-folded recipe calls); '
              'no sign of contraction')
        return 0
    if per > a.max_per_sqrt:
        print(f'UNEXPLAINED: {rest} FMAs left, {per:.1f} per sqrt > {a.max_per_sqrt}: possible contraction, '
              'stop and tell the user (SKILL.md Phase 6)')
        return 1
    print(f'OK: remainder {per:.1f} FMAs per sqrt, no sign of contraction')
    return 0


if __name__ == '__main__':
    sys.exit(main())
