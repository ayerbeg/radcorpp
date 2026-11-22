#!/usr/bin/env python3
"""Read a radcor --details file and print, per spectrum, the RC factor range and the mean tail fractions.
usage: python3 read_details.py details.txt"""
import sys
import numpy as np

names = None
rows = []
for line in open(sys.argv[1]):
    if line.startswith("#spec"):
        names = line[1:].split()          # column names
    elif line.strip() and not line.startswith("#"):
        rows.append([float(x) for x in line.split()])
d = dict(zip(names, np.array(rows).T))

for s in np.unique(d["spec"]).astype(int):
    m = d["spec"] == s
    rc = d["RC"][m]
    tails = (d["tail_Es"][m] + d["tail_Ep"][m]) / d["xs_in"][m]
    print(f"E0={d['E0'][m][0]:7.1f} MeV  points={m.sum():4d}  "
          f"RC=[{rc.min():.3f},{rc.max():.3f}]  mean tail fraction={tails.mean():.3f}")
