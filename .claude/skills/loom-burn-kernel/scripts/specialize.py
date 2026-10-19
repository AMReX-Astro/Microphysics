#!/usr/bin/env python3
"""Phase 3 starting point: C++ Insights' specialized source for each reachable function.

usage: uv run --no-project specialize.py --inv DIR [--test-dir unit_test/test_react]
           [--insights insights] [--llvm DIR]

--llvm is the LLVM whose clang resource dir (builtin headers) insights is
given; it must match the LLVM insights links against. Default: $LLVM_BIN if
set, else the LLVM that the insights binary links (otool -L / ldd; Ubuntu's
/lib/x86_64-linux-gnu/libLLVM-18.so.1 maps to /usr/lib/llvm-18/bin), else
the LLVM major version from `insights --version`, else the dir of `clang` on
PATH, else /opt/homebrew/opt/llvm/bin (Homebrew's insights links llvm@21).
--test-dir: relative paths that don't exist from the cwd are taken relative
to the Microphysics root (the repo holding this skill).

Installing insights: macOS `brew install cppinsights`. Ubuntu has no
package: build https://github.com/andreasfertig/cppinsights with CMake at
the tag matching the distro LLVM. On Ubuntu 24.04 (tested): tag v_20.1, apt
clang-20 llvm-20-dev libclang-20-dev, and
-DINSIGHTS_LLVM_CONFIG=/usr/lib/llvm-20/bin/llvm-config (SKILL.md Setup).
`main` needs LLVM >= 21; v_17.0 against Ubuntu's LLVM 18 fails to link and,
patched, crashes on main.cpp. Or use a GitHub release binary and pass --llvm for the LLVM version it was built with. The
LLVM must match: a mismatched -resource-dir gives builtin-header errors.

Needs an inventory directory (inventory.py --out DIR): it reuses flags.txt and
the reachable `definition` rows of hits.tsv. Runs `insights` on main.cpp and
writes
  DIR/specialized.cpp          the whole TU as Insights prints it
  DIR/specialized/<name>.cpp   every concrete definition of <name>: template
                               instantiations (template<> ... f<Args>(...)) and
                               non-template definitions; uninstantiated
                               template patterns are skipped
  DIR/specialized/INDEX.md     name -> source location, instantiation count

Insights resolves templates, `if constexpr`, `auto`, implicit conversions,
user-defined literals and operator calls into explicit C++. It does NOT
make types visible at call sites (`a.operator()(i, j)` does not say whether
`a` is a 1-based Array2D or a 0-based Array4), so index rewrites still need
the declarations.

LITERALS ON x86-64 (e.g. an Ubuntu VM): Insights prints literals from the
host's long double. On arm64 that is binary64, so
`amrex::operator""_rt(0.10000000000000001L)` is exactly the double the source
meant. On x86-64 long double is 80-bit, so a printed literal may carry digits
that round differently: always copy constants from the original source, never
from this output.

Lambdas are printed as closure classes, one per instantiation, so
`constexpr_for` loops (the generic network RHS/Jacobian in networks/rhs.H,
dgefa/dgesl) come out fully unrolled: for aprox13 `jac` is ~93k lines. Read
loop structure from the original source and use this output for the resolved
types, branches and constants. (Fully unrolled LU lost 7-13% in the
chemistry campaigns; don't translate the unrolled form by default.)

Use this output to read and translate from, never compile it as the kernel.
"""
import argparse
import collections
import os
import re
import shutil
import subprocess
import sys


DEFAULT_LLVM = '/opt/homebrew/opt/llvm/bin'
SKILL = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.abspath(os.path.join(SKILL, '..', '..', '..'))


def insights_llvm(insights):
    """bin dir of the LLVM the insights binary links against, or None."""
    exe = shutil.which(insights)
    if not exe:
        return None
    tool = ['otool', '-L'] if sys.platform == 'darwin' else ['ldd']
    try:
        out = subprocess.run([*tool, exe], capture_output=True, text=True).stdout
    except OSError:
        return None
    m = re.search(r'(\S+)/lib/libLLVM[^/\s]*\.(?:dylib|so)', out)
    if m and os.path.isfile(os.path.join(m.group(1), 'bin', 'clang')):
        return os.path.join(m.group(1), 'bin')
    # Debian/Ubuntu: /lib/x86_64-linux-gnu/libLLVM-18.so.1 or libLLVM.so.18.1
    m = re.search(r'libLLVM(?:-(\d+)|\.so\.(\d+))', out)
    if m:
        d = f'/usr/lib/llvm-{m.group(1) or m.group(2)}/bin'
        if os.path.isfile(os.path.join(d, 'clang')):
            return d
    # statically linked: ask insights which LLVM it was built with
    try:
        ver = subprocess.run([exe, '--version'], capture_output=True, text=True).stdout
    except OSError:
        ver = ''
    m = re.search(r'LLVM version (\d+)', ver)
    if m:
        d = f'/usr/lib/llvm-{m.group(1)}/bin'
        if os.path.isfile(os.path.join(d, 'clang')):
            return d
    return None


def clang_on_path():
    exe = shutil.which('clang')
    return os.path.dirname(os.path.realpath(exe)) if exe else None


def run_insights(a, test, flags):
    sdk = subprocess.run(['xcrun', '--show-sdk-path'], capture_output=True,
                         text=True).stdout.strip() if sys.platform == 'darwin' else ''
    try:
        rd = subprocess.run([f'{a.llvm}/clang', '-print-resource-dir'], capture_output=True,
                            text=True).stdout.strip()
    except OSError:
        sys.exit(f'specialize: no clang in {a.llvm}; pass --llvm (or set LLVM_BIN) to the '
                 'LLVM bin dir that insights was built against, e.g. /usr/lib/llvm-20/bin')
    extra = ['-resource-dir', rd] if rd else []
    if sdk:  # Homebrew's insights does not find the SDK or libc++ on its own
        extra += ['-isysroot', sdk, '-isystem', f'{sdk}/usr/include/c++/v1']
    try:
        r = subprocess.run([a.insights, 'main.cpp', '--', *extra, *flags], cwd=test,
                           capture_output=True, text=True)
    except OSError:
        sys.exit(f'specialize: cannot run {a.insights!r} (C++ Insights); see the install notes '
                 'at the top of this script, or pass --insights PATH')
    if r.returncode:
        sys.exit('specialize: insights failed:\n' + '\n'.join(
            l for l in r.stderr.splitlines() if 'error' in l)[:4000])
    return r.stdout


def body_end(lines, start):
    """Index of the line closing the brace block that opens at or after `start`."""
    depth, opened = 0, False
    for i in range(start, len(lines)):
        t = re.sub(r'"(\\.|[^"\\])*"|\'(\\.|[^\'\\])*\'', '""', lines[i])
        t = re.sub(r'//.*', '', t)
        for ch in t:
            if ch == '{':
                depth += 1
                opened = True
            elif ch == '}':
                depth -= 1
        if opened and depth <= 0:
            return i
        if not opened and t.rstrip().endswith(';'):
            return None  # a declaration, not a definition
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--inv', required=True)
    ap.add_argument('--test-dir', default='unit_test/test_react')
    ap.add_argument('--insights', default='insights')
    ap.add_argument('--llvm', default=None,
                    help='LLVM bin dir matching insights (default: $LLVM_BIN, else detected '
                         'from the insights binary, else clang on PATH, else '
                         '/opt/homebrew/opt/llvm/bin)')
    a = ap.parse_args()
    a.llvm = (a.llvm or os.environ.get('LLVM_BIN') or insights_llvm(a.insights)
              or clang_on_path() or DEFAULT_LLVM)
    inv = os.path.abspath(a.inv)
    test = a.test_dir
    if not (os.path.isabs(test) or os.path.isdir(test)):
        test = os.path.join(ROOT, test)
    test = os.path.abspath(test)
    if not os.path.isfile(os.path.join(test, 'main.cpp')):
        sys.exit(f'specialize: no main.cpp in --test-dir {test}')
    flags = open(os.path.join(inv, 'flags.txt')).read().split()
    defs = collections.defaultdict(set)
    for row in open(os.path.join(inv, 'hits.tsv')).read().splitlines()[1:]:
        cat, loc, fn, reach, _ = (row.split('\t') + [''])[:5]
        if cat == 'definition' and reach == 'Y' and re.match(r'^~?\w+$', fn):
            defs[fn].add(loc)
    if not defs:
        sys.exit('specialize: no reachable definitions in hits.tsv; rerun inventory.py')

    text = run_insights(a, test, flags)
    open(os.path.join(inv, 'specialized.cpp'), 'w').write(text)
    lines = text.split('\n')
    head = re.compile(r'^\s*[^=;(){}]*?\b(' + '|'.join(map(re.escape, sorted(defs)))
                      + r')\s*(<.*?>)?\s*\(')
    found = collections.defaultdict(list)
    i = 0
    attrs = re.compile(r'__attribute__\s*\(\(.*?\)\)\s*')
    while i < len(lines):
        # Insights prints every definition as a one-line header, then "{" alone
        m = head.match(attrs.sub('', lines[i]))
        if (not m or i + 1 >= len(lines) or lines[i + 1].strip() != '{'
                or re.match(r'\s*(return|if|for|while|else|switch)\b', lines[i])):
            i += 1
            continue
        prev = next((lines[j].strip() for j in range(i - 1, max(i - 4, -1), -1)
                     if lines[j].strip()), '')
        pattern = prev.startswith('template') and prev != 'template<>'
        end = body_end(lines, i)
        if end is None:
            i += 1
            continue
        if not pattern:
            first = i
            while first > 0 and (lines[first - 1].startswith(('template<>', '#ifdef INSIGHTS',
                                                              '/* First instantiated'))):
                first -= 1
            found[m.group(1)].append('\n'.join(lines[first:end + 1]))
        i = end + 1

    outdir = os.path.join(inv, 'specialized')
    os.makedirs(outdir, exist_ok=True)
    idx = ['# Specialized definitions (C++ Insights)', '',
           '| function | source | definitions | lines |', '|---|---|---|---|']
    total = 0
    for fn in sorted(defs):
        blocks = found.get(fn, [])
        n = sum(b.count('\n') + 1 for b in blocks)
        total += n
        if blocks:
            open(os.path.join(outdir, f'{fn}.cpp'), 'w').write('\n\n'.join(blocks) + '\n')
        idx.append(f'| {fn} | {", ".join(sorted(defs[fn]))} | {len(blocks)} | {n} |')
    missing = [fn for fn in defs if not found.get(fn)]
    idx += ['', f'{total} lines in {len(defs) - len(missing)} files.',
            f'Not found ({len(missing)}): ' + ', '.join(sorted(missing)), '']
    open(os.path.join(outdir, 'INDEX.md'), 'w').write('\n'.join(idx))
    print(f'specialize: {len(defs) - len(missing)}/{len(defs)} functions, {total} lines -> {outdir}')


if __name__ == '__main__':
    main()
