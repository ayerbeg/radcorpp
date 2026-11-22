#!/usr/bin/env python3
"""Unfold the E94-010 spectra from numpy arrays with the Python interface (python/radcor.py).
usage: python3 examples/python_unfold.py      (from the suite top directory, after building)"""
import csv
import os
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "python"))
import radcor  # noqa: E402

# Get the E94-010 points as arrays (here by converting the shipped ele-input to CSV).
with tempfile.TemporaryDirectory() as tmp:
    table = os.path.join(tmp, "e94010.csv")
    subprocess.run([radcor.find_executable(), "convert", "unpol", "--input",
                    os.path.join(HERE, "..", "data", "e94010", "unpol", "ele-input"), "--to", table, "-q"], check=True)
    rows = list(csv.DictReader(l for l in open(table) if not l.startswith("#")))

spectra = {}
for r in rows:
    s = spectra.setdefault(r["E0"], dict(E0=float(r["E0"]), nu=[], xs=[], stat=[], syst=[],
                                         t_before=float(r["t_before"]), t_after=float(r["t_after"])))
    for k in ("nu", "xs", "stat", "syst"):
        s[k].append(float(r[k]))
spectra = [dict(s, **{k: np.array(s[k]) for k in ("nu", "xs", "stat", "syst")}) for s in spectra.values()]

res = radcor.unfold(spectra, Z=2, A=3, theta=15.5)  # clean settings, double precision
for e0 in np.unique(res["E0"]):
    m = res["E0"] == e0
    print(f"E0={e0:7.1f} MeV  points={m.sum():4d}  RC=[{res['rc'][m].min():.3f}, {res['rc'][m].max():.3f}]")
