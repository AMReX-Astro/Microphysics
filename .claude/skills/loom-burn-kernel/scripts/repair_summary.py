#!/usr/bin/env python3
"""Summarize loom-compile repair traces (LOOM_REPAIR lines).

usage: repair_summary.py [TRACE ...] [--json]
       repair_summary.py -h

TRACE is a file with the stderr of a loom-compile built with
patches/trace-repair.patch and run with LOOM_TRACE_REPAIR=1 (other lines are
ignored); with no TRACE, or '-', stdin is read. Each LOOM_REPAIR line is
  LOOM_REPAIR round= limit= values= schedule_errors= allocation_errors=
              spill_plans= rematerialized_operands= split_operands= round_seconds=
(the short spellings remat= and split= are accepted too). round=0 starts a
new spill batch; spill_plans=0 means the repair converged.

Prints: total rounds, number of batches, per batch its starting spill_plans,
rounds and seconds (and values/limit at the start), total repair seconds,
first and last values, final rematerialized/split operand counts, and whether
the last round converged. --json prints the same as one JSON object.
Python 3 standard library only; runs locally or on the GPU host.
"""
import argparse
import json
import re
import sys

ALIASES = {'remat': 'rematerialized_operands', 'split': 'split_operands'}
PAIR = re.compile(r'(\w+)=(\S+)')


def parse(lines):
    rounds = []
    for line in lines:
        i = line.find('LOOM_REPAIR')
        if i < 0:
            continue
        r = {}
        for k, v in PAIR.findall(line[i + len('LOOM_REPAIR'):]):
            k = ALIASES.get(k, k)
            try:
                r[k] = float(v) if k == 'round_seconds' or '.' in v else int(v)
            except ValueError:
                r[k] = v
        if 'round' in r:
            rounds.append(r)
    return rounds


def summarize(rounds):
    batches = []
    for r in rounds:
        if r['round'] == 0 or not batches:
            batches.append([])
        batches[-1].append(r)
    first, last = (rounds[0], rounds[-1]) if rounds else ({}, {})
    return {
        'rounds': len(rounds),
        'batches': len(batches),
        'batch_table': [{
            'batch': i + 1,
            'start_spill_plans': b[0].get('spill_plans'),
            'end_spill_plans': b[-1].get('spill_plans'),
            'start_values': b[0].get('values'),
            'limit': b[0].get('limit'),
            'rounds': len(b),
            'seconds': round(sum(r.get('round_seconds', 0.0) for r in b), 1),
        } for i, b in enumerate(batches)],
        'repair_seconds': round(sum(r.get('round_seconds', 0.0) for r in rounds), 1),
        'first_values': first.get('values'),
        'last_values': last.get('values'),
        'final_rematerialized_operands': last.get('rematerialized_operands'),
        'final_split_operands': last.get('split_operands'),
        'final_schedule_errors': last.get('schedule_errors'),
        'final_allocation_errors': last.get('allocation_errors'),
        'converged': bool(rounds) and last.get('spill_plans') == 0,
    }


def show(s, name):
    print(f'== {name}')
    if not s['rounds']:
        print('no LOOM_REPAIR lines (built with patches/trace-repair.patch and run '
              'with LOOM_TRACE_REPAIR=1?)')
        return
    print(f"rounds={s['rounds']} batches={s['batches']} repair_seconds={s['repair_seconds']} "
          f"converged={'yes' if s['converged'] else 'NO'}")
    print(f"values first={s['first_values']} last={s['last_values']}  "
          f"final remat={s['final_rematerialized_operands']} split={s['final_split_operands']} "
          f"schedule_errors={s['final_schedule_errors']} "
          f"allocation_errors={s['final_allocation_errors']}")
    cols = ['batch', 'start_spill_plans', 'end_spill_plans', 'start_values', 'limit',
            'rounds', 'seconds']
    rows = [[str(b[c]) for c in cols] for b in s['batch_table']]
    w = [max(len(c), *(len(r[i]) for r in rows)) for i, c in enumerate(cols)]
    print('  '.join(c.rjust(w[i]) for i, c in enumerate(cols)))
    for r in rows:
        print('  '.join(v.rjust(w[i]) for i, v in enumerate(r)))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0],
                                 epilog=__doc__[__doc__.index('TRACE is'):],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('trace', nargs='*', default=['-'], help="trace file(s); '-' = stdin")
    ap.add_argument('--json', action='store_true', help='print JSON instead of a table')
    a = ap.parse_args()
    out = {}
    for t in a.trace:
        if t == '-':
            s = summarize(parse(sys.stdin))
        else:
            with open(t, errors='replace') as f:
                s = summarize(parse(f))
        out['stdin' if t == '-' else t] = s
    if a.json:
        print(json.dumps(out if len(out) > 1 else next(iter(out.values())), indent=2))
    else:
        for name, s in out.items():
            show(s, name)
    return 0 if all(s['rounds'] for s in out.values()) else 1


if __name__ == '__main__':
    sys.exit(main())
