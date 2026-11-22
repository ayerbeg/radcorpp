# radcorpp — radiative corrections for inclusive electron scattering

`radcor` is a C++17 translation of **nradcor**, the Fortran code that removes radiative effects from
measured inclusive electron-scattering spectra. It takes measured cross sections
d²σ/dΩdE′ at one scattering angle and several beam energies, and returns the Born (non-radiative)
cross section at the same points. It can also do the reverse: radiate a Born model ("folding").

The method is the iterative unfolding of Mo and Tsai (1969), Sec. IV, with the radiating formula of
Stein *et al.* (1975), Eq. A82. It covers internal bremsstrahlung, external bremsstrahlung in the
target and windows, and ionisation loss.

**The physics and the algorithm are not ours.** They come from the authors credited
[below](#credits). This repository is a translation, documentation and validation of their work.

---

## What this translation adds

- **Same results as the Fortran.** In `--precision float` mode, with the default settings, the output is byte-identical
  to the Fortran on the shipped E94-010 inputs, both unpolarized and polarized. This is checked by
  `ctest` on every build.
- **Every known defect of the Fortran is documented, measured and switchable.** By default the C++
  reproduces each one, so results stay comparable with the Fortran; `--set legacy=0` turns them all off,
  and each can also be switched separately.
  For new work use `--precision double --set legacy=0`.
- **Transparent output.** For every point it writes the two radiative-tail integrals, the soft+virtual
  factor, the radiated model, the Born result, RC = σ_meas/σ_Born and the errors. These go to a text
  table (`--details`) and/or a ROOT tree (`--root`), together with the full run settings.
- **Flexible input.** It reads the original Fortran files, a CSV table with named columns, a ROOT tree,
  or numpy arrays from Python. `radcor convert` translates between these formats without loss, and
  `radcor check` validates an input before it is used.
- **A manual** ([`doc/radcor_manual.pdf`](doc/radcor_manual.pdf), 70 pages) covering:
  - the physics, with equations cited from the papers;
  - what goes in and what comes out;
  - every option;
  - every known defect, with its measured size;
  - the validation.

## Quick start

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release     # -DRADCOR_WITH_ROOT=OFF to build without ROOT
cmake --build build -j8
ctest --test-dir build                              # parity with the Fortran, input formats, Python

# reproduce the Fortran on the E94-010 input (byte-identical to data/e94010/unpol/expected/)
build/radcor unpol --dir data/e94010/unpol --out e94010.dat

# recommended settings for new data
build/radcor unpol --input mydata.csv --precision double --set legacy=0 \
     --set collisional_loss=0 --set finite_acc=0 --details details.txt --root details.root

build/radcor check   unpol --input mydata.csv                   # validate an input
build/radcor convert unpol --input ele-input --to mydata.csv    # Fortran file -> CSV (or .root)
build/radcor -h                                                 # all options
```

From Python (no compiled bindings needed; it runs the executable):

```python
import sys; sys.path.insert(0, "python"); import radcor
spectra = [dict(E0=862.0,  nu=nu1, xs=xs1, stat=st1, t_before=7.184e-3, t_after=4.166e-2),
           dict(E0=1716.9, nu=nu2, xs=xs2, stat=st2, t_before=7.922e-3, t_after=7.492e-2)]
r = radcor.unfold(spectra, Z=2, A=3, theta=15.5)     # double precision, legacy bugs off
r["xs_born"], r["rc"], r["tail_Es"], r["stat"]       # numpy arrays, one entry per point
```

**Requirements:**
- a C++17 compiler (tested with g++ 16.2) and CMake ≥ 3.16;
- optional: ROOT (Tree, RIO; tested with 6.40), Python 3 with numpy, and gfortran to rebuild the original Fortran.

## The input, in one table

Whatever the format, a run needs the following:

| Level | Quantities |
|---|---|
| run | target Z and A, scattering angle θ [deg], number of passes |
| spectrum | beam energy E₀ (spectra in increasing E₀), material before and after the vertex [radiation lengths], normalisation, RC uncertainty fraction |
| point | ν = E₀ − E′ [MeV] (increasing), radiated cross section (any unit), statistical and systematic errors |

**Before you use it on new data:**
- Subtract the elastic (and, for nuclei, quasi-elastic) radiative tail first. The code does not compute it.
- The spectra should start close to the quasi-elastic or elastic region. Chapter 3 of the manual explains why.

## Layout

```
apps/radcor.cpp            command-line front end (run, check, convert)
include/radcor/Engine.h    the unfolding engine: radiating formula, Born model, iteration
include/radcor/Input.h     the input, independent of format; formats; validation
include/radcor/Run.h       drivers: input -> engine -> outputs in the Fortran formats
include/radcor/Options.h   every switch and parameter, with its default
include/radcor/Interp.h    spline and polynomial interpolation of the output stage
src/                       formats, drivers, option parsing, ROOT input/output
python/radcor.py           Python interface
fortran/                   the original Fortran sources (see Credits), Makefile
data/e94010/               E94-010 inputs and the Fortran reference outputs
examples/                  reading the outputs; using the library from C++ and Python
tests/                     ctest: Fortran parity, input formats, Python, Fortran rebuild
doc/                       the manual; K. Slifer's description of his modifications
```

## Known limitations of the method

The manual (chapter "Known issues") gives the full list, with measured sizes. In short:

- **No elastic or quasi-elastic radiative tail.** It must be subtracted from the input beforehand.
- **Peaking approximation, Stein form.** It differs from Miller's form by 10–20 % at ν ≳ E_s/3. Slifer
  noted this in the source.
- **Model edges.** The Born model is interpolated at fixed ν/E. It is zero below the first point of each
  spectrum, constant above the last, and scaled below the lowest beam energy.
- **Not included:** Coulomb corrections and two-photon exchange.
- **E94-010 specifics.** The collisional-loss shift and the finite-acceptance factor are specific to that
  experiment. They are switches (`collisional_loss`, `finite_acc`) and should be off for other data.
- **Polarized version.** It applies the unpolarized kernels to cross-section differences.

## Credits

The radiative-correction method, the original program, and its modifications are the work of:

- **L. W. Mo and Y. S. Tsai.** The theory and the iterative unfolding method:
  *Radiative corrections to elastic and inelastic ep and μp scattering*,
  Rev. Mod. Phys. **41**, 205 (1969).
- **S. Stein *et al.*** (SLAC). The radiating formula used by the code (Eq. A82):
  *Electron scattering at 4° with energies of 4.5–20 GeV*, Phys. Rev. D **12**, 1884 (1975).
- **G. Miller.** The basis of the original program: *Inelastic Electron Scattering at Large Angles*,
  Ph.D. thesis, SLAC-R-129 (1971), Eq. A.V.5.
- **L. Whitney (SLAC).** Author of the original unfolding program `radcor`. The version distributed with
  nradcor is an adaptation made at Saclay ("PROGRAMME WHITNEY ADAPTE AU SUN"); the person who made it
  is not named in the source. It is kept unchanged in
  [`fortran/whitney/`](fortran/whitney).
- **K. Slifer (Temple University, 2003–2005).** Author of **nradcor**, the modified program translated
  here, written for JLab experiment E94-010:
  - cleaned code and a new input format;
  - the Stein et al. formalism for the soft factor and the integrands;
  - the collisional-loss shift and the finite-acceptance correction;
  - the polarized version.

  His sources are in [`fortran/unpolarized/`](fortran/unpolarized) and
  [`fortran/polarized/`](fortran/polarized). His description of the changes is in
  [`doc/slifer_newcodes.pdf`](doc/slifer_newcodes.pdf).
- **Seonho Choi.** The interpolation routine `INTERP` ("Interpol_1dim") used in the output stage.
- **W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P. Flannery.** The cubic-spline routines
  `SPLINE`/`SPLINT`, from *Numerical Recipes* (Cambridge University Press).
- **The JLab E94-010 collaboration** (Hall A). The ³He spectra and finite-acceptance tables in
  [`data/e94010/`](data/e94010).
  - The acceptance tables were generated with **QFS**, the quasi-free scattering code of
    J. W. Lightbody Jr. and J. S. O'Connell.
  - The polarized example comes from a Hall C kinematics (5.775 GeV, 13.07°). Its first seven spectra are
    model input from MAID2003 plus a Hall C model.

**The translation:** C. Ayerbe (JLab). The C++ translation, its validation and the manual were written with the help of
Claude (Anthropic).

If you use this code, please cite Mo and Tsai and Stein *et al.*, and acknowledge Whitney's and Slifer's
programs.

## Licence

No licence has been chosen yet for the C++ code. The Fortran programs in `fortran/` remain the work of
their authors and are included unchanged for reference and validation.
