#!/usr/bin/env python3
"""Phase 2 inventory: what the burn reaches, and which Loom-hostile constructs it uses.

usage (from the Microphysics root or anywhere):
  uv run --no-project inventory.py --out DIR [--test-dir unit_test/test_react]
      [--root do_react] [--llvm DIR] -- MAKEVAR=VALUE ...
  --llvm: LLVM bin dir with clang++/clang-query (default $LLVM_BIN, else the
  dir of clang-query on PATH, else the newest /usr/lib/llvm-*/bin with
  clang-query (Ubuntu: apt install clang-N clang-tools-N), else Homebrew's
  /opt/homebrew/opt/llvm/bin; specialize.py instead follows insights' LLVM)
  --test-dir: relative paths that don't exist from the cwd are taken relative
  to the Microphysics root (the repo holding this skill).

Prerequisites: AMReX checked out and AMREX_HOME exported (make -n fails
otherwise; its output is kept in DIR/make-n.log), GNU make and g++.

Steps (all local, CPU only, no Loom tools):
  1. Generate the configuration's headers (network_properties.H,
     extern_parameters.H, AMReX_Config.H) into DIR/build via
     TMP_BUILD_DIR, so the test's own tmp_build_dir is never touched or
     reused with another network's headers. main.cpp itself is only
     syntax-checked.
  2. Take main.cpp's -D/-I/-std flags from `make -n`.
  3. clang's debug.DumpCallGraph + closure from --root -> reachable.txt.
  4. clang-query matchers over Microphysics sources -> hits.tsv, then
     summary.md (counts for functions reachable from --root).

Caveats: reachability is by printed name, so overloads and template
specializations sharing a name are merged (an upper bound). The construct
list mirrors references/loom-frontend-rules.md and must be kept in sync with
the Phase 1 probes; constructs not matched here (ternary nesting depth, base
classes, operator-overloaded types other than autodiff) still need a look.
"""
import argparse
import collections
import glob
import os
import re
import shlex
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True   # no __pycache__ in the (possibly read-only) skill dir
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import reachable  # noqa: E402

SKILL = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def _microphysics_root():
    """$MICROPHYSICS_HOME, else the repo this skill sits in
    (<repo>/.claude/skills/<skill> or <repo>/.agents/skills/<skill>), else the
    git checkout of the current directory (a user-level skill install)."""
    def ok(d):
        return d and os.path.isdir(os.path.join(d, 'unit_test', 'test_react'))
    env = os.environ.get('MICROPHYSICS_HOME')
    if ok(env):
        return os.path.abspath(env)
    here = os.path.abspath(os.path.join(SKILL, '..', '..', '..'))
    if ok(here):
        return here
    try:
        top = subprocess.run(['git', 'rev-parse', '--show-toplevel'], capture_output=True,
                             text=True).stdout.strip()
    except OSError:
        top = ''
    if ok(top):
        return top
    sys.exit('cannot find the Microphysics tree: set MICROPHYSICS_HOME or run from inside it')


ROOT = _microphysics_root()

SRC = (r'(^|/)(networks|integration|EOS|rates|screening|neutrinos|interfaces|'
       r'constants|util|nse_solver|nse_tabular)/|react_zones\\.H$')
# vendored libraries are counted at their call sites (dual-var, math:*), not inside
F = (f'isExpansionInFileMatching("{SRC}"), '
     'unless(isExpansionInFileMatching("/util/(autodiff|gcem)/"))')
FN = ('hasAncestor(functionDecl(unless(isImplicit()), unless(cxxMethodDecl('
      'ofClass(cxxRecordDecl(isLambda()))))).bind("fn"))')
MATH = ['exp', 'log', 'pow', 'log10', 'log2', 'exp2', 'exp10', 'sqrt', 'cbrt',
        'tanh', 'sinh', 'cosh', 'sin', 'cos', 'tan', 'atan', 'atan2', 'asinh',
        'expm1', 'log1p', 'erf', 'erfc', 'fmod', 'floor', 'ceil', 'round',
        'trunc', 'fabs', 'abs', 'powi', 'fast_exp', 'fast_log', 'fast_pow']

MATCHERS = [
    ('ref-param', f'parmVarDecl({F}, hasType(referenceType()), {FN})'),
    ('ref-local', f'varDecl({F}, hasLocalStorage(), unless(parmVarDecl()), '
                  f'hasType(referenceType()), {FN})'),
    ('ref-return', f'functionDecl({F}, isDefinition(), unless(isImplicit()), '
                   'unless(cxxMethodDecl(ofClass(cxxRecordDecl(isLambda())))), '
                   'returns(referenceType()))'
                   '.bind("fn")'),
    ('member-call', f'cxxMemberCallExpr({F}, {FN}, '
                    'callee(cxxMethodDecl().bind("callee")))'),
    ('operator-call', f'cxxOperatorCallExpr({F}, {FN}, '
                      'callee(functionDecl().bind("callee")))'),
    ('lambda', f'lambdaExpr({F}, {FN})'),
    ('switch', f'switchStmt({F}, {FN})'),
    ('return-in-loop', f'returnStmt({F}, {FN}, hasAncestor(stmt(anyOf('
                       'forStmt(), whileStmt(), doStmt(), cxxForRangeStmt()))))'),
    ('range-for', f'cxxForRangeStmt({F}, {FN})'),
    ('for-ever', f'forStmt({F}, {FN}, unless(hasCondition(expr())))'),
    ('uninit-local', f'varDecl({F}, {FN}, hasLocalStorage(), unless(parmVarDecl()), '
                     'unless(hasInitializer(anything())))'),
    ('global-var', f'declRefExpr({F}, {FN}, to(varDecl(hasGlobalStorage(), '
                   'unless(isConstexpr()), unless(hasType(isConstQualified())), '
                   'optionally(hasDeclContext(namespaceDecl().bind("ns"))))'
                   '.bind("var")))'),
    ('constexpr-array', f'declRefExpr({F}, {FN}, to(varDecl(hasGlobalStorage(), '
                        'isConstexpr(), hasType(hasCanonicalType(arrayType())), '
                        'optionally(hasDeclContext(namespaceDecl().bind("ns"))))'
                        '.bind("var")))'),
    ('dual-var', f'varDecl({F}, {FN}, hasType(hasUnqualifiedDesugaredType('
                 'recordType(hasDeclaration(cxxRecordDecl('
                 'matchesName("autodiff::detail::Dual")))))))'),
    ('definition', f'functionDecl({F}, isDefinition(), unless(isImplicit()), '
                   'unless(cxxMethodDecl(ofClass(cxxRecordDecl(isLambda())))))'
                   '.bind("fn")'),
] + [(f'math:{m}', f'callExpr({F}, {FN}, callee(functionDecl(hasName("{m}"))'
                   '.bind("callee")))') for m in MATH]

DIAG = re.compile(r'^(.+?):(\d+):(\d+): note: "(\w+)" binds here')
SKIP_CALLS = {'if', 'for', 'while', 'switch', 'return', 'sizeof', 'alignof',
              'decltype', 'alignas', '__attribute__', 'noexcept', 'static_assert',
              'requires'}


def run(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


def default_llvm():
    if os.environ.get('LLVM_BIN'):
        return os.environ['LLVM_BIN']
    cq = shutil.which('clang-query')
    if cq:
        return os.path.dirname(os.path.realpath(cq))
    def ver(d):
        m = re.search(r'llvm-(\d+)', d)
        return int(m.group(1)) if m else 0
    cands = sorted((d for d in glob.glob('/usr/lib/llvm-*/bin')
                    if os.path.exists(os.path.join(d, 'clang-query'))), key=ver)
    if cands:
        return cands[-1]
    return '/opt/homebrew/opt/llvm/bin'


def resolve_test_dir(d):
    if os.path.isabs(d) or os.path.isdir(d):
        return os.path.abspath(d)
    return os.path.join(ROOT, d)


def strip_angles(t):
    out, depth = '', 0
    for ch in t:
        if ch == '<':
            depth += 1
        elif ch == '>' and depth:
            depth -= 1
        elif depth == 0:
            out += ch
    return out


class Source:
    def __init__(self, cwd):
        self.cwd, self.cache = cwd, {}

    def text_from(self, path, line, col, nlines=12):
        p = path if os.path.isabs(path) else os.path.join(self.cwd, path)
        if p not in self.cache:
            try:
                self.cache[p] = open(p, errors='replace').read().split('\n')
            except OSError:
                self.cache[p] = []
        lines = self.cache[p][line - 1:line - 1 + nlines]
        if not lines:
            return ''
        lines[0] = lines[0][col - 1:]
        return '\n'.join(lines)

    def function_name(self, path, line, col):
        t = self.text_from(path, line, col, 20)
        t = re.sub(r'//[^\n]*', '', t)
        while True:  # drop template headers, balancing <>
            m = re.search(r'\btemplate\s*<', t)
            if not m:
                break
            i, depth = m.end(), 1
            while i < len(t) and depth:
                depth += {'<': 1, '>': -1}.get(t[i], 0)
                i += 1
            t = t[:m.start()] + t[i:]
        for m in re.finditer(r'(operator\s*(?:\(\)|[^\s(]+)|~?[A-Za-z_]\w*)\s*\(', t):
            if m.group(1) not in SKIP_CALLS:
                return re.sub(r'\s+', '', m.group(1))
        return '?'

    def var_name(self, path, line, col):
        t = strip_angles(self.text_from(path, line, col, 3))
        t = re.split(r'[\[=;{(,]', t, maxsplit=1)[0]
        ids = re.findall(r'[A-Za-z_]\w*', t)
        return ids[-1] if ids else '?'

    def ns_name(self, path, line, col):
        m = re.match(r'\s*(?:inline\s+)?namespace\s+([\w:]+)', self.text_from(path, line, col, 2))
        return m.group(1) if m else ''


def short(name):
    """Last unqualified component of a call-graph name, lambdas mapped to their owner."""
    name = re.sub(r'\(.*', '', name)          # "f(args)::(lambda)..." -> "f"
    depth, out = 0, ''
    for ch in name:                            # drop template arguments
        depth += {'<': 1, '>': -1}.get(ch, 0) if not out.endswith('operator') else 0
        if depth == 0 and ch not in '<>':
            out += ch
    return out.split('::')[-1]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--out', required=True)
    ap.add_argument('--test-dir', default='unit_test/test_react')
    ap.add_argument('--root', action='append')
    ap.add_argument('--llvm', default=None,
                    help='LLVM bin dir with clang++ and clang-query (default: $LLVM_BIN, '
                         'else clang-query on PATH, else newest /usr/lib/llvm-*/bin, '
                         'else /opt/homebrew/opt/llvm/bin)')
    ap.add_argument('--hip-arch', default='gfx90a',
                    help="also take the USE_HIP build's -D set (via a stub hipconfig/hipcc "
                         "if no ROCm is installed) and report how it differs; '' to skip")
    ap.add_argument('makevars', nargs='*')
    a = ap.parse_args()
    roots = a.root or ['do_react']
    a.llvm = a.llvm or default_llvm()
    missing = [t for t in ('clang++', 'clang-query') if not os.access(os.path.join(a.llvm, t), os.X_OK)]
    if missing:
        sys.exit(f'inventory: {" and ".join(missing)} not found in {a.llvm}; install LLVM '
                 '(Ubuntu: apt install clang-N clang-tools-N, then --llvm /usr/lib/llvm-N/bin; '
                 'macOS: brew install llvm) or pass --llvm / set LLVM_BIN')
    out = os.path.abspath(a.out)
    os.makedirs(out, exist_ok=True)
    test = resolve_test_dir(a.test_dir)
    if not os.path.isfile(os.path.join(test, 'main.cpp')):
        sys.exit(f'inventory: no main.cpp in --test-dir {test}')
    build = os.path.join(out, 'build')
    mk = ['make', f'TMP_BUILD_DIR={build}'] + a.makevars

    # 1-2. compile line, then generate headers by building main.o syntax-only
    dryrun = run(mk + ['-n'], cwd=test)
    dry = dryrun.stdout
    open(os.path.join(out, 'make-n.log'), 'w').write(dryrun.stdout + dryrun.stderr)
    line = next((l for l in dry.splitlines() if re.search(r' -c main\.cpp( |$)', l)), None)
    if line is None:
        sys.exit(f'inventory: no main.cpp compile line in `make -n` output (see {out}/make-n.log; '
                 'is AMREX_HOME exported and pointing at an AMReX checkout?)')
    toks = shlex.split(line)
    obj = toks[toks.index('-o') + 1]
    gen = run(mk + ['XTRA_CXXFLAGS=-fsyntax-only', obj], cwd=test)
    open(os.path.join(out, 'make.log'), 'w').write(gen.stdout + gen.stderr)
    if not any('network_properties.H' in f for _, _, fs in os.walk(build) for f in fs):
        sys.exit(f'inventory: header generation failed, see {out}/make.log')
    flags, i = [], 1
    while i < len(toks):
        t = toks[i]
        if t.startswith(('-D', '-U', '-I', '-std=')):
            flags.append(t)
        elif t in ('-include', '-isystem', '-iquote'):
            flags += toks[i:i + 2]
            i += 1
        i += 1
    open(os.path.join(out, 'flags.txt'), 'w').write('\n'.join(flags) + '\n')

    hipdiff = hip_define_diff(a, test, build, out, flags) if a.hip_arch else None

    # 3. call graph and closure
    cg = run([f'{a.llvm}/clang++', *flags, '--analyze', '-Xanalyzer',
              '-analyzer-checker=debug.DumpCallGraph', '-o', '/dev/null', 'main.cpp'],
             cwd=test)
    cgpath = os.path.join(out, 'callgraph.txt')
    open(cgpath, 'w').write(cg.stdout + cg.stderr)
    r = run([sys.executable, reachable.__file__, cgpath, *roots])
    if r.returncode:
        sys.exit(r.stderr.strip())
    reach = r.stdout.split('\n')
    collapsed = sorted({re.sub(r'std::integral_constant<int, \d+>', 'N', n)
                        for n in reach if n})
    open(os.path.join(out, 'reachable.txt'), 'w').write('\n'.join(collapsed) + '\n')
    reach_short = {short(n) for n in reach if n}

    # 4. matchers
    qfile = os.path.join(out, 'queries.txt')
    with open(qfile, 'w') as q:
        q.write('set traversal AsIs\nset output diag\nset print-matcher true\n')
        for _, m in MATCHERS:
            q.write(f'match {m}\n')
    cq = run([f'{a.llvm}/clang-query', '-f', qfile, 'main.cpp', '--', *flags], cwd=test)
    open(os.path.join(out, 'clang-query.out'), 'w').write(cq.stdout + cq.stderr)

    hits = relabel(cq.stdout, Source(test), reach_short)

    with open(os.path.join(out, 'hits.tsv'), 'w') as h:
        h.write('category\tlocation\tfunction\treachable\tdetail\n')
        for row in hits:
            h.write('\t'.join(map(str, row)) + '\n')
    write_summary(out, hits, collapsed, roots, flags, a.makevars, hipdiff)
    print(f'inventory: {len(collapsed)} reachable names, {len(hits)} distinct hits -> {out}/summary.md')


def hip_define_diff(a, test, build, out, cpu_flags):
    """-D/-U set of the USE_HIP build's main.cpp compile line, diffed against the CPU one.

    AMReX's hip.mak needs `hipconfig --path` and `hipcc --version` at parse time.
    Without ROCm, a stub directory provides both; `make -n` never runs hipcc.
    """
    env = dict(os.environ)
    if not shutil.which('hipconfig'):
        stub = os.path.join(out, 'hipstub')
        os.makedirs(os.path.join(stub, 'bin'), exist_ok=True)
        for name, body in [('hipconfig', f'case "$1" in --path) echo {stub};; '
                                         '--compiler) echo clang;; *) echo 6.4.0;; esac'),
                           ('hipcc', 'echo "HIP version: 6.4.0-stub"')]:
            path = os.path.join(stub, 'bin', name)
            open(path, 'w').write(f'#!/bin/sh\n{body}\n')
            os.chmod(path, 0o755)
        env['PATH'] = os.path.join(stub, 'bin') + os.pathsep + env['PATH']
    r = subprocess.run(['make', '-n', f'TMP_BUILD_DIR={build}-hip', 'USE_HIP=TRUE',
                        f'AMD_ARCH={a.hip_arch}'] + a.makevars,
                       cwd=test, env=env, capture_output=True, text=True)
    open(os.path.join(out, 'make-n-hip.log'), 'w').write(r.stdout + r.stderr)
    line = next((l for l in r.stdout.splitlines() if re.search(r' -c main\.cpp( |$)', l)), None)
    if line is None:
        return {'error': f'no HIP compile line; see {out}/make-n-hip.log'}
    hip = {t for t in shlex.split(line) if t.startswith(('-D', '-U'))}
    cpu = {t for t in cpu_flags if t.startswith(('-D', '-U'))}
    fp = [t for t in shlex.split(line) if t.startswith(('-ffp-contract', '-ffast-math', '-fno-fast'))]
    return {'cpu_only': sorted(cpu - hip), 'hip_only': sorted(hip - cpu), 'hip_fp': fp}


def relabel(text, src, reach_short):
    """Parse clang-query diag output; blocks are introduced by 'Matcher:' lines."""
    hits, seen = [], set()
    names = [n for n, _ in MATCHERS]
    idx, binds = -1, {}

    def flush():
        if idx < 0 or 'root' not in binds:
            return
        cat = names[idx]
        fpath, ln, col = binds['root']
        fn = src.function_name(*binds['fn']) if 'fn' in binds else '?'
        detail = ''
        if 'var' in binds:
            ns = src.ns_name(*binds['ns']) if 'ns' in binds else ''
            detail = (ns + '::' if ns else '') + src.var_name(*binds['var'])
        elif 'callee' in binds:
            cf = binds['callee'][0]
            kind = ('autodiff' if 'autodiff' in cf else
                    'array' if re.search(r'Array', cf) else
                    'std' if re.search(r'/c\+\+/|/usr/|SDKs', cf) else '')
            detail = (kind + ':' if kind else '') + src.function_name(*binds['callee'])
        key = (cat, fpath, ln, col, detail)
        if key not in seen:
            seen.add(key)
            hits.append((cat, f'{os.path.normpath(fpath)}:{ln}', fn,
                         'Y' if fn in reach_short else 'N', detail))

    printed_matcher = False
    for l in text.splitlines():
        s = l.strip()
        if s.startswith('Matcher:'):
            # print-matcher echoes the matcher before the results and again
            # before the count; only the first echo starts a new block
            if not printed_matcher:
                flush()
                binds = {}
                idx += 1
                printed_matcher = True
            continue
        if re.match(r'^\d+ match(es)?\.$', s) or s == 'No matches.':
            flush()
            binds = {}
            printed_matcher = False
            continue
        if l.startswith('Match #'):
            flush()
            binds = {}
            continue
        m = DIAG.match(l)
        if m:
            binds[m.group(4)] = (m.group(1), int(m.group(2)), int(m.group(3)))
    flush()
    return hits


def write_summary(out, hits, reach, roots, flags, makevars, hipdiff=None):
    R = [h for h in hits if h[3] == 'Y']
    by_cat = collections.Counter(h[0] for h in R)
    by_cat_all = collections.Counter(h[0] for h in hits)
    s = ['# Loom inventory', '',
         f'make vars: `{" ".join(makevars)}`  ', f'roots: {", ".join(roots)}  ',
         f'reachable names (collapsed): {len(reach)} — see reachable.txt  ',
         'defines: `' + ' '.join(f for f in flags if f.startswith('-D')) + '`', '',
         ] + hip_section(hipdiff) + [
         '## Construct hits (distinct source locations)', '',
         '| category | in reachable functions | all Microphysics |', '|---|---|---|']
    for c, _ in MATCHERS:
        if by_cat_all[c] and c != 'definition':
            s.append(f'| {c} | {by_cat[c]} | {by_cat_all[c]} |')
    defs = sorted({(h[2], h[1]) for h in R if h[0] == 'definition'},
                  key=lambda d: d[1])
    s += ['', f'## Reachable Microphysics function definitions ({len(defs)})', '']
    s += [f'- `{fn}` {loc}' for fn, loc in defs]
    s += ['', '## Globals read on reachable paths (become BurnParams / BurnTables)', '']
    g = collections.Counter(h[4] for h in R if h[0] == 'global-var')
    s += [f'- `{v}` ({n})' for v, n in sorted(g.items())] or ['- none']
    s += ['', '## constexpr arrays read on reachable paths', '']
    # AMReX's terminal colour/font tables are read only by Print/Abort paths
    noise = ('FGColor::', 'BGColor::', 'Font::', 'amrex::ResetDisplay')
    g = collections.Counter(h[4] for h in R if h[0] == 'constexpr-array'
                            and not h[4].startswith(noise))
    s += [f'- `{v}` ({n})' for v, n in sorted(g.items())] or ['- none']
    s += ['', '## Math calls on reachable paths (callee kind:name)', '']
    g = collections.Counter(f'{h[0][5:]} -> {h[4]}' for h in R if h[0].startswith('math:'))
    s += [f'- {v}: {n}' for v, n in sorted(g.items())] or ['- none']
    s += ['', '## Per-function hits (reachable only)', '',
          '| function | ' + ' | '.join(c for c, _ in MATCHERS if by_cat[c] and not c.startswith('math:') and c != 'definition') + ' |',
          '|---|' + '---|' * sum(1 for c, _ in MATCHERS if by_cat[c] and not c.startswith('math:') and c != 'definition')]
    cats = [c for c, _ in MATCHERS if by_cat[c] and not c.startswith('math:') and c != 'definition']
    per = collections.defaultdict(collections.Counter)
    for h in R:
        per[h[2]][h[0]] += 1
    for fn in sorted(per, key=lambda f: -sum(per[f][c] for c in cats)):
        if any(per[fn][c] for c in cats):
            s.append(f'| {fn} | ' + ' | '.join(str(per[fn][c] or '') for c in cats) + ' |')
    s += ['', 'These results use the -D set of the build above. A USE_HIP build differs',
          '(e.g. HIP_CONSTEXPR tables, no ALLOW_JACOBIAN_CACHING, AMREX_USE_GPU branches);',
          'diff its `make -n` defines against these before trusting the counts.',
          '', 'Reachability is by unqualified name (upper bound). Locations, functions and',
          'details for every hit are in hits.tsv; raw output in clang-query.out.', '']
    open(os.path.join(out, 'summary.md'), 'w').write('\n'.join(s))


def hip_section(d):
    if d is None:
        return []
    if 'error' in d:
        return ['## CPU vs HIP defines', '', d['error'], '']
    s = ['## CPU vs HIP defines', '',
         'The analysis below uses the CPU set. Defines that differ change what the burn',
         'does on GPU; the CPU reference for Phase 4 must match the HIP build on the',
         'physics ones (e.g. strip -DALLOW_JACOBIAN_CACHING).', '',
         '- CPU only: ' + (' '.join(f'`{x}`' for x in d['cpu_only']) or 'none'),
         '- HIP only: ' + (' '.join(f'`{x}`' for x in d['hip_only']) or 'none'),
         '- HIP fp flags on the compile line: ' + (' '.join(f'`{x}`' for x in d['hip_fp'])
                                                  or 'none (hipcc default contracts: pass XTRA_CXXFLAGS=-ffp-contract=off)'),
         '']
    return s


if __name__ == '__main__':
    main()
