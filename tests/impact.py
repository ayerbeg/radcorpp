#!/usr/bin/env python3
"""impact.py A B : per spectrum, relative change (A/B-1) of XSf (col 4): mean, median|.|, max|.| and nu at max."""
import sys, statistics as st, collections
def load(f):
    d = {}
    for l in open(f):
        if l.startswith('#') or not l.strip(): continue
        v = l.split(); d[(v[0], float(v[1]))] = float(v[3])
    return d
a, b = load(sys.argv[1]), load(sys.argv[2])
by = collections.defaultdict(list)
for k in a:
    if k in b and b[k] != 0: by[k[0]].append((k[1], a[k] / b[k] - 1))
for e, r in by.items():
    rr = [x for _, x in r]; nu, mx = max(r, key=lambda t: abs(t[1]))
    print(f"  {e:>7}: mean {st.mean(rr):+.2e}  median|.| {st.median(map(abs, rr)):.1e}  max {mx:+.2e} at nu={nu:g}")
