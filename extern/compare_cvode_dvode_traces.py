#!/usr/bin/env python3
import math
import re
import sys
from collections import defaultdict

KV = re.compile(r"([A-Za-z0-9_]+)=\s*([-+]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?|[-+]?\.\d+(?:[Ee][+-]?\d+)?)")


def parse(path, tags):
    out = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            sline = line.lstrip()
            tag = None
            for t in tags:
                if sline.startswith(t):
                    tag = t
                    break
            if tag is None:
                continue
            body = sline[len(tag):].strip()
            if not body:
                continue
            event = body.split()[0]
            # Keep control-channel events distinct so they do not
            # inflate counts for core integration events.
            if tag.endswith("_CTL] "):
                event = f"CTL:{event}"
            kvs = {k: float(v) for k, v in KV.findall(sline)}
            out[event].append(kvs)
    return out


def summarize(lhs, rhs, name_lhs, name_rhs):
    report = []
    all_events = sorted(set(lhs) | set(rhs))
    for ev in all_events:
        n1 = len(lhs.get(ev, []))
        n2 = len(rhs.get(ev, []))
        n = min(n1, n2)
        keys = set()
        for i in range(n):
            keys.update(lhs[ev][i].keys())
            keys.update(rhs[ev][i].keys())
        key_stats = []
        for k in sorted(keys):
            diffs = []
            for i in range(n):
                if k in lhs[ev][i] and k in rhs[ev][i]:
                    diffs.append(abs(lhs[ev][i][k] - rhs[ev][i][k]))
            if diffs:
                key_stats.append((k, max(diffs), diffs[0]))
        key_stats.sort(key=lambda x: x[1], reverse=True)
        report.append((ev, n1, n2, n, key_stats[:8]))
    return report


def canonical_count(events, event, required_keys):
    return sum(1 for kv in events.get(event, []) if all(k in kv for k in required_keys))


def main():
    if len(sys.argv) != 3:
        print("usage: compare_cvode_dvode_traces.py <cvode_log> <dvode_log>")
        return 2

    cv = parse(sys.argv[1], ("[CVODE] ", "[CVODE_CTL] "))
    dv = parse(sys.argv[2], ("[DVODE] ", "[DVODE_CTL] "))

    rep = summarize(cv, dv, "CVODE", "DVODE")

    for ev, ncv, ndv, n, top in rep:
        print(f"EVENT {ev}: CVODE={ncv} DVODE={ndv} compared={n}")
        for k, m, f in top:
            print(f"  {k}: max_abs_diff={m:.6e} first_abs_diff={f:.6e}")

    # quick overall line for core events used in both traces
    core = ["PRE", "POST", "REJECT", "ACCEPT", "ORDER_APPLY", "ORDER_DECIDE", "dvnlsd:"]
    print("SUMMARY core_event_counts:")
    for ev in core:
        print(f"  {ev}: CVODE={len(cv.get(ev, []))} DVODE={len(dv.get(ev, []))}")

    # Canonical counts filter duplicate log variants that share the same event token.
    # Example: PRE appears in multiple lines; only the line with tn/H is the step PRE event.
    canonical_keys = {
        "PRE": ("tn", "H"),
        "POST": ("ACNRM", "DSM"),
        "REJECT": ("DSM", "ETA"),
        "ACCEPT": ("t", "hu", "nq"),
        "ORDER_APPLY": ("ETA", "NEWQ"),
        "ORDER_DECIDE": ("DSM",),
        # Count one nonlinear-iteration event per iteration using the DCON line.
        "dvnlsd:": ("DCON",),
    }
    print("SUMMARY canonical_core_event_counts:")
    for ev in core:
        keys = canonical_keys.get(ev, ())
        ncv = canonical_count(cv, ev, keys)
        ndv = canonical_count(dv, ev, keys)
        print(f"  {ev}: CVODE={ncv} DVODE={ndv}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
