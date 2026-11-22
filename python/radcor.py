"""Python interface to radcor.

Runs the radcor executable on spectra given as Python/numpy arrays and returns the per-point
results (the --details table) as numpy arrays. No compiled bindings are needed: the spectra are
written to a CSV input file (format "csv" of the manual), radcor is run, and its details table
is read back.

    import radcor
    spectra = [dict(E0=862.0, nu=nu1, xs=xs1, stat=st1, t_before=7.184e-3, t_after=4.166e-2),
               dict(E0=1716.9, nu=nu2, xs=xs2, t_before=7.922e-3, t_after=7.492e-2), ...]
    r = radcor.unfold(spectra, Z=2, A=3, theta=15.5)
    r["xs_born"], r["rc"], r["tail_Es"], ...      # one entry per point, all spectra

The executable is found from $RADCOR_EXE, then ../build/radcor next to this file, then PATH.
"""
import os
import shutil
import subprocess
import tempfile

try:
    import numpy as np
except ImportError:  # numpy is optional: plain lists are returned without it
    np = None

# Spectrum-level keys and their defaults (as in the CSV format); point-level keys are arrays.
_SPECTRUM_KEYS = {"t_target": 0.0, "t_before": 0.0, "t_after": 0.0, "rc_sys": 0.0, "norm": 1.0,
                  "loss_s": None, "loss_p": None, "acc_table": "", "title": ""}
_POINT_KEYS = ("nu", "xs", "stat", "syst")

# Recommended settings for new work: all known Fortran bugs off, no E94-010 specifics.
CLEAN = {"legacy": 0, "collisional_loss": 0, "finite_acc": 0}


def find_executable():
    """Path of the radcor executable ($RADCOR_EXE, ../build/radcor, or PATH)."""
    exe = os.environ.get("RADCOR_EXE")
    if exe:
        return exe
    here = os.path.dirname(os.path.abspath(__file__))
    local = os.path.join(here, "..", "build", "radcor")
    if os.access(local, os.X_OK):
        return local
    exe = shutil.which("radcor")
    if exe:
        return exe
    raise FileNotFoundError("radcor executable not found; set RADCOR_EXE")


def _fmt(v):
    """Number as text that reads back exactly (empty for None)."""
    return "" if v is None else repr(float(v))


def write_csv(path, spectra, Z, A, theta, kind="unpol", nitt=10):
    """Write spectra to a radcor CSV input file.

    spectra: list of dicts, one per beam energy. Required keys: E0, nu, xs (arrays for nu, xs).
    Optional: stat, syst (arrays); t_target, t_before, t_after [radiation lengths], rc_sys, norm,
    loss_s, loss_p [MeV], acc_table, title. Polarized (kind="pol"): xs, stat are pairs
    (xs=(xs1, xs2), stat=(st1, st2)) on the same nu.
    """
    cols = ["E0", "config", "nu", "xs", "stat", "syst"] + list(_SPECTRUM_KEYS)
    with open(path, "w") as f:
        f.write(f"#! kind={kind} nitt={int(nitt)} Z={_fmt(Z)} A={_fmt(A)} theta={_fmt(theta)}\n")
        f.write(",".join(cols) + "\n")
        for s in spectra:
            nu = list(s["nu"])
            nconf = 2 if kind == "pol" else 1
            for k in range(nconf):
                xs = s["xs"][k] if nconf == 2 else s["xs"]
                st = s.get("stat")
                st = (st[k] if nconf == 2 else st) if st is not None else [0.0] * len(nu)
                sy = s.get("syst") if nconf == 1 else None
                sy = sy if sy is not None else [0.0] * len(nu)
                spec = [s.get(key, d) for key, d in _SPECTRUM_KEYS.items()]
                spec = [str(v).replace(",", " ") if isinstance(v, str) else _fmt(v) for v in spec]
                for j in range(len(nu)):
                    row = [_fmt(s["E0"]), str(k + 1), _fmt(nu[j]), _fmt(xs[j]), _fmt(st[j]), _fmt(sy[j])] + spec
                    f.write(",".join(row) + "\n")


def read_details(path):
    """Read a radcor --details table into a dict of columns (numpy arrays if numpy is present).
    The run settings (the '#!' header lines) are returned under the key 'settings'."""
    names, rows, settings = None, [], []
    with open(path) as f:
        for line in f:
            if line.startswith("#!"):
                settings.append(line[2:].rstrip())
            elif line.startswith("#spec"):
                names = line[1:].split()
            elif line.strip() and not line.startswith("#"):
                rows.append([float(x) for x in line.split()])
    if names is None:
        raise ValueError(f"{path}: no column header")
    cols = list(zip(*rows)) if rows else [[] for _ in names]
    out = {n: (np.array(c) if np else list(c)) for n, c in zip(names, cols)}
    out["rc"] = out["RC"]  # lower-case alias, as the ROOT branch
    out["settings"] = "\n".join(settings)
    return out


def run(kind, input, options=None, precision="double", details=None, out=None, root=None,
        data=None, fmt=None, exe=None, quiet=True):
    """Run radcor on an existing input file (any format) and return the details table.

    options: dict of --set keys (e.g. {"legacy": 0, "delta": 5}). out/root: optional output files.
    """
    exe = exe or find_executable()
    with tempfile.TemporaryDirectory() as tmp:
        det = details or os.path.join(tmp, "details.txt")
        cmd = [exe, kind, "--input", os.path.abspath(input), "--precision", precision,
               "--details", os.path.abspath(det), "--out", os.path.abspath(out or os.path.join(tmp, "out.dat"))]
        if data:
            cmd += ["--data", os.path.abspath(data)]
        if fmt:
            cmd += ["--format", fmt]
        if root:
            cmd += ["--root", os.path.abspath(root)]
        for k, v in (options or {}).items():
            cmd += ["--set", f"{k}={v}"]
        if quiet:
            cmd.append("-q")
        p = subprocess.run(cmd, capture_output=True, text=True)
        if p.returncode != 0:
            raise RuntimeError(f"radcor failed ({p.returncode}):\n{p.stderr}")
        return read_details(det)


def unfold(spectra, Z, A, theta, nitt=10, options=None, precision="double", kind="unpol", **kw):
    """Unfold spectra given as arrays (see write_csv); returns the details table as a dict.
    Default options: CLEAN (all legacy bugs off, no E94-010 specifics); pass options to change."""
    opts = dict(CLEAN)
    opts.update(options or {})
    with tempfile.TemporaryDirectory() as tmp:
        csv = os.path.join(tmp, "input.csv")
        write_csv(csv, spectra, Z, A, theta, kind=kind, nitt=nitt)
        return run(kind, csv, options=opts, precision=precision, **kw)


def radiate(spectra, Z, A, theta, options=None, precision="double", **kw):
    """Radiate Born spectra (folding): returns the details table; the radiated model is 'rad_model'.
    Same input as unfold(), with the Born cross section in xs."""
    opts = dict(CLEAN)
    opts.update({"radiate": 1, "nitt": 1, "radiate_overwrite": 0})
    opts.update(options or {})
    return unfold(spectra, Z, A, theta, nitt=1, options=opts, precision=precision, **kw)
