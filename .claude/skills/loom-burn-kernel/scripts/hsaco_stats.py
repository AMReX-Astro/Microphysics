#!/usr/bin/env python3
"""Per-kernel resources and arguments of an AMDGPU code object (.hsaco).

usage: hsaco_stats.py HSACO [--llvm DIR] [--report REPORT.json] [--json]
       hsaco_stats.py --notes-text NOTES.txt [--report REPORT.json] [--json]
       hsaco_stats.py --report REPORT.json
       hsaco_stats.py -h

For every kernel in the code-object metadata (from `llvm-readelf --notes`)
prints: name, vgpr/agpr/sgpr counts, private_segment_fixed_size (scratch per
work-item), group_segment_fixed_size (LDS), kernarg_segment_size, wavefront
size, and the kernel arguments (name, offset, size, value_kind) in order.
The metadata is parsed with a small built-in YAML-subset reader (no PyYAML).

llvm-readelf is looked up in --llvm DIR, then $LLVM, then ROCm's
($ROCM_PATH/llvm/bin, /opt/rocm/llvm/bin, /opt/rocm*/lib/llvm/bin), then PATH
last: an older distro llvm-readelf may not decode newer AMDGPU code-object
metadata. --notes-text parses a saved `llvm-readelf --notes` output instead
of running it.

--report REPORT.json (written by loom-compile --compile-report=summary
--compile-report-output=FILE; .json or .json.gz) additionally prints every
scalar whose flattened dotted path contains spill, scratch, register,
pressure, time, second, occupancy, private, allocation, repair or limiting.
Per-row tables (lists such as schedule_band_summary_rows: 58,483 rows for a
full burn, an 84 MB report) are summarized as row counts unless --all is
given. Gzip the report after reading it.

Python 3 standard library only; meant to run on the GPU host as well.
"""
import argparse
import glob
import gzip
import json
import os
import re
import shutil
import subprocess
import sys


# ------------------------------------------------------------ YAML subset

def _scalar(v):
    v = v.strip()
    if not v:
        return None
    if v[0] in '\'"' and v[-1] == v[0] and len(v) >= 2:
        return v[1:-1]
    if v[0] == '[' and v[-1] == ']':
        return [_scalar(x) for x in v[1:-1].split(',') if x.strip()]
    if re.fullmatch(r'[-+]?\d+', v):
        return int(v)
    if re.fullmatch(r'0x[0-9a-fA-F]+', v):
        return int(v, 16)
    if v in ('true', 'false'):
        return v == 'true'
    return v


def _lines(text):
    out = []
    for raw in text.splitlines():
        if not raw.strip() or raw.strip().startswith('#'):
            continue
        out.append([len(raw) - len(raw.lstrip(' ')), raw.strip()])
    return out


def _parse(ls, i, ind):
    """Parse the block starting at ls[i] with indentation ind -> (obj, next i)."""
    if ls[i][1].startswith('- ') or ls[i][1] == '-':
        seq = []
        while i < len(ls) and ls[i][0] == ind and (ls[i][1].startswith('- ') or ls[i][1] == '-'):
            body = ls[i][1][2:] if ls[i][1] != '-' else ''
            if not body:
                if i + 1 < len(ls) and ls[i + 1][0] > ind:
                    obj, i = _parse(ls, i + 1, ls[i + 1][0])
                else:
                    obj, i = None, i + 1
                seq.append(obj)
            elif re.match(r'^[^\'"\[]+?:(\s|$)', body):
                # "- key: value" starts a mapping indented by the dash width
                ls[i] = [ind + 2, body]
                obj, i = _parse(ls, i, ind + 2)
                seq.append(obj)
            else:
                seq.append(_scalar(body))
                i += 1
        return seq, i
    m = {}
    while i < len(ls) and ls[i][0] == ind and not ls[i][1].startswith('- '):
        mm = re.match(r'^(.+?):(?:\s+(.*))?$', ls[i][1])
        if not mm:
            i += 1
            continue
        key, val = mm.group(1).strip(), (mm.group(2) or '').strip()
        i += 1
        if val:
            m[key] = _scalar(val)
        elif i < len(ls) and (ls[i][0] > ind or (ls[i][0] == ind and ls[i][1].startswith('- '))):
            m[key], i = _parse(ls, i, ls[i][0])
        else:
            m[key] = None
    return m, i


def parse_yaml_docs(text):
    """All '---' ... '...' documents in readelf --notes output."""
    docs, cur = [], None
    for raw in text.splitlines():
        s = raw.strip()
        if s == '---':
            cur = []
        elif s == '...' and cur is not None:
            ls = _lines('\n'.join(cur))
            if ls:
                base = min(l[0] for l in ls)
                for l in ls:
                    l[0] -= base
                docs.append(_parse(ls, 0, 0)[0])
            cur = None
        elif cur is not None:
            cur.append(raw)
    return docs


# ------------------------------------------------------------ code object

def find_readelf(llvm):
    cands = []
    for d in (llvm, os.environ.get('LLVM')):
        if d:
            cands.append(os.path.join(d, 'llvm-readelf'))
    rocm = os.environ.get('ROCM_PATH')
    if rocm:
        cands.append(os.path.join(rocm, 'llvm', 'bin', 'llvm-readelf'))
    cands += ['/opt/rocm/llvm/bin/llvm-readelf'] + sorted(
        glob.glob('/opt/rocm*/lib/llvm/bin/llvm-readelf') + glob.glob('/opt/rocm/*/lib/llvm/bin/llvm-readelf'))
    w = shutil.which('llvm-readelf')
    if w:
        cands.append(w)
    for c in cands:
        if os.access(c, os.X_OK):
            return c
    sys.exit('hsaco_stats: llvm-readelf not found (use --llvm DIR or set LLVM)')


def kernels(docs):
    out = []
    for d in docs:
        if isinstance(d, dict):
            for k in d.get('amdhsa.kernels') or []:
                if isinstance(k, dict):
                    k = {key.lstrip('.'): v for key, v in k.items()}
                    k['args'] = [{x.lstrip('.'): v for x, v in a.items()} if isinstance(a, dict) else a
                                 for a in (k.get('args') or [])]
                    out.append(k)
    return out


FIELDS = [('vgpr', 'vgpr_count'), ('agpr', 'agpr_count'), ('sgpr', 'sgpr_count'),
          ('scratch', 'private_segment_fixed_size'), ('lds', 'group_segment_fixed_size'),
          ('kernarg', 'kernarg_segment_size'), ('wave', 'wavefront_size')]


def show_kernels(ks):
    if not ks:
        print('no amdhsa.kernels in the code-object metadata')
    for k in ks:
        print(f"kernel {k.get('name')}  (symbol {k.get('symbol')})")
        print('  ' + '  '.join(f"{lbl}={k.get(key)}" for lbl, key in FIELDS))
        extra = {x: k[x] for x in ('sgpr_spill_count', 'vgpr_spill_count', 'max_flat_workgroup_size',
                                   'uses_dynamic_stack') if x in k}
        if extra:
            print('  ' + '  '.join(f'{x}={v}' for x, v in extra.items()))
        args = k.get('args') or []
        print(f'  args ({len(args)}):')
        for a in args:
            a = a if isinstance(a, dict) else {}
            print(f"    offset={a.get('offset')!s:>4} size={a.get('size')!s:>3} "
                  f"{a.get('value_kind')!s:<26} {a.get('name', '')}"
                  + (f" ({a['address_space']})" if 'address_space' in a else ''))


# ------------------------------------------------------------ compile report

REPORT_WORDS = ('spill', 'scratch', 'register', 'pressure', 'time', 'second',
                'occupancy', 'private', 'allocation', 'repair', 'limiting')


def flatten(o, path=''):
    if isinstance(o, dict):
        for k, v in o.items():
            yield from flatten(v, f'{path}.{k}' if path else str(k))
    elif isinstance(o, list):
        for i, v in enumerate(o):
            yield from flatten(v, f'{path}[{i}]')
    else:
        yield path, o


def report_items(path, all_rows=False):
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as f:
        rep = json.load(f)
    items = []
    tables = {}
    for p, v in flatten(rep):
        if '[' in p and not all_rows:
            # one entry of a per-row table: count the rows, don't print them
            head = p[:p.index('[')]
            tables.setdefault(head, set()).add(p[:p.index(']') + 1])
            continue
        if any(w in p.lower() for w in REPORT_WORDS):
            items.append((p, v))
    items += [(f'{head} (rows)', len(rows)) for head, rows in sorted(tables.items())
              if len(rows) > 50]
    return items


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0],
                                 epilog=__doc__[__doc__.index('For every kernel'):],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('hsaco', nargs='?', help='code object (.hsaco)')
    ap.add_argument('--llvm', help='dir with llvm-readelf')
    ap.add_argument('--notes-text', help='saved `llvm-readelf --notes` output to parse instead')
    ap.add_argument('--report', help='loom-compile JSON compile report (.json or .json.gz)')
    ap.add_argument('--all', action='store_true', help='with --report: also print every per-row table entry')
    ap.add_argument('--json', action='store_true', help='print JSON')
    a = ap.parse_args()
    if not (a.hsaco or a.notes_text or a.report):
        ap.print_help()
        return 2
    result = {}
    if a.hsaco or a.notes_text:
        if a.notes_text:
            with open(a.notes_text, errors='replace') as f:
                text = f.read()
        else:
            r = subprocess.run([find_readelf(a.llvm), '--notes', a.hsaco],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               universal_newlines=True)  # python >= 3.6
            if r.returncode:
                sys.exit(f'hsaco_stats: llvm-readelf failed: {r.stderr.strip()}')
            text = r.stdout
        result['kernels'] = kernels(parse_yaml_docs(text))
    if a.report:
        result['report'] = dict(report_items(a.report, a.all))
    if a.json:
        print(json.dumps(result, indent=2, default=str))
        return 0
    if 'kernels' in result:
        show_kernels(result['kernels'])
    if 'report' in result:
        print(f'compile report {a.report}:')
        if not result['report']:
            print('  (no spill/scratch/register/pressure/time keys)')
        for p, v in result['report'].items():
            print(f'  {p} = {v}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
