#!/usr/bin/env python3
"""Compare two nradcor-format result files column by column (matching on the first key columns).
usage: compare.py A B [nkey]   -> relative differences per spectrum/column"""
import sys, collections
def load(f, nkey):
    d = collections.OrderedDict()
    for l in open(f):
        if l.startswith('#') or not l.strip(): continue
        v = l.split(); d[tuple(v[:nkey])] = [float(x) for x in v[nkey:]]
    return d
nkey = int(sys.argv[3]) if len(sys.argv) > 3 else 2
a, b = load(sys.argv[1], nkey), load(sys.argv[2], nkey)
print(f"rows: {len(a)} vs {len(b)}; only in A: {len(set(a)-set(b))}; only in B: {len(set(b)-set(a))}")
by = collections.defaultdict(list)
for k in a:
    if k in b: by[k[0]].append((a[k], b[k]))
worst = 0.0
for e, rows in by.items():
    ncol = len(rows[0][0]); out = []
    for c in range(ncol):
        m = 0.0; nd = 0
        for x, y in rows:
            den = max(abs(x[c]), abs(y[c]))
            if den > 0:
                r = abs(x[c] - y[c]) / den; m = max(m, r); nd += r > 0
        out.append(f"{m:.1e}({nd})"); worst = max(worst, m)
    print(f"{e:>8} n={len(rows):4d}  max rel diff per column (#differing): " + " ".join(out))
print(f"WORST {worst:.2e}")
