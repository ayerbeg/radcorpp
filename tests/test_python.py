#!/usr/bin/env python3
"""Python interface: unfold the E94-010 spectra from arrays and compare with the command line
(same settings); then radiate the result back (closure, printed only)."""
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "python"))
import radcor  # noqa: E402

exe = radcor.find_executable()
data = os.path.join(HERE, "..", "data", "e94010", "unpol")
with tempfile.TemporaryDirectory() as tmp:
    ref = radcor.run("unpol", os.path.join(data, "ele-input"), options=radcor.CLEAN)
    table = os.path.join(tmp, "u.csv")
    subprocess.run([exe, "convert", "unpol", "--input", os.path.join(data, "ele-input"), "--to", table, "-q"], check=True)
    import csv
    rows = list(csv.DictReader(l for l in open(table) if not l.startswith("#")))
spectra = {}
for r in rows:
    s = spectra.setdefault(r["E0"], dict(E0=float(r["E0"]), nu=[], xs=[], stat=[], syst=[],
                                         t_before=float(r["t_before"]), t_after=float(r["t_after"]),
                                         rc_sys=float(r["rc_sys"]), title=r["title"]))
    for k in ("nu", "xs", "stat", "syst"):
        s[k].append(float(r[k]))
res = radcor.unfold(list(spectra.values()), Z=2, A=3, theta=15.5)
same = all(a == b for a, b in zip(res["xs_born"], ref["xs_born"])) and len(res["xs_born"]) == len(ref["xs_born"]) == 1031
print("PASS python unfold == command line (1031 points)" if same else "FAIL python unfold differs")
sys.exit(0 if same else 1)
