#!/usr/bin/env python3
"""Reachable-function closure over clang's debug.DumpCallGraph output.

usage: reachable.py CALLGRAPH.txt ROOT [ROOT ...] > reachable.txt

Clang prints one "Function: <name> calls: <callee> <callee> ..." record per
function, names without parameter types, and wraps long records over several
lines. Names can contain spaces ("unsigned long", "Dual<double, X>") and
angle brackets in operator names, so callees are split on spaces outside
<...> and (...) after masking operator spellings.

Template specializations and overloads that share a printed name are merged,
so the closure over-approximates; treat it as an upper bound.
"""
import argparse
import collections
import re
import sys

OPS = re.compile(r'operator(<=>|<<=|>>=|<<|>>|<=|>=|->\*|->|<|>|\(\))')


def records(path):
    cur = None
    for line in open(path, errors='replace'):
        line = line.rstrip('\n')
        if line.lstrip().startswith('Function: '):
            if cur is not None:
                yield cur
            cur = line.strip()
        elif cur is not None:
            cur += ' ' + line.strip()
    if cur is not None:
        yield cur


def split_names(s):
    masked = []
    s = OPS.sub(lambda m: masked.append(m.group(0)) or '\0', s)
    out, cur, depth = [], '', 0
    for ch in s:
        if ch in '<(':
            depth += 1
        elif ch in '>)':
            depth -= 1
        if ch == ' ' and depth == 0:
            if cur:
                out.append(cur)
            cur = ''
        else:
            cur += ch
    if cur:
        out.append(cur)
    it = iter(masked)
    return [re.sub('\0', lambda _: next(it), n) for n in out]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0],
                                 epilog=__doc__[__doc__.index('Clang prints'):],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('callgraph', help="clang debug.DumpCallGraph output (CALLGRAPH.txt)")
    ap.add_argument('roots', nargs='+', metavar='ROOT', help='function name(s) as printed in the call graph')
    a = ap.parse_args()
    path, roots = a.callgraph, a.roots
    graph = collections.defaultdict(set)
    for rec in records(path):
        m = re.match(r'Function: (.*?) calls:(.*)$', rec)
        if m:
            graph[m.group(1).strip()] |= set(split_names(m.group(2).strip()))
    missing = [r for r in roots if r not in graph]
    if missing:
        sys.exit('reachable.py: root(s) not in call graph: ' + ', '.join(missing))
    seen, todo = set(), list(roots)
    while todo:
        f = todo.pop()
        if f not in seen:
            seen.add(f)
            todo.extend(graph.get(f, ()))
    for f in sorted(seen):
        print(f)


if __name__ == '__main__':
    main()
