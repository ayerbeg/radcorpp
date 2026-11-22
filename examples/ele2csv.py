#!/usr/bin/env python3
"""Convert an ele-input file to a flat CSV table, one row per point:
    E0,nu,xs,stat,syst,t_target,t_before,t_after
(xs and stat multiplied by RNORME). Used by engine_from_csv.cpp to show that the engine
does not depend on the ele-input format.
usage: python3 ele2csv.py ele-input > points.csv"""
import sys
L = open(sys.argv[1]).read().splitlines()
num = lambda s: s.split("!")[0].replace(",", " ").split()
ne = int(num(L[1])[0])
thick = [list(map(float, num(L[4 + i]))) for i in range(ne)]   # TX TIW TFW RC_SYS RNORME
print("E0,nu,xs,stat,syst,t_target,t_before,t_after")
pos = 4 + ne
for i in range(ne):
    e0, _, _, n = num(L[pos + 1]); n = int(float(n))
    tx, tiw, tfw, _, rn = thick[i]
    for l in L[pos + 2: pos + 2 + n]:
        nu, xs, st, sy = map(float, num(l)[:4])
        print(f"{e0},{nu},{xs * rn:.6e},{st * rn:.6e},{sy:.6e},{tx},{tiw},{tfw}")
    pos += 2 + n
