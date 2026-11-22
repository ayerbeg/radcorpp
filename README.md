# radcor

C++17 translation of the radiative-correction unfolding code `nradcor` (L. Whitney, SLAC;
K. Slifer, JLab E94-010). Physics: Mo & Tsai, Rev. Mod. Phys. 41, 205 (1969), Sec. IV, with the
radiating formula of Stein et al., Phys. Rev. D 12, 1884 (1975), Eq. A82.

By default (`--precision float`, legacy switches on) the output is byte-identical to the Fortran.
For new work use `--precision double --set legacy=0`. The manual (doc/radcor_manual.pdf) documents
the physics, the input and output files, every option and every known defect.

## Build and test
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release      # -DRADCOR_WITH_ROOT=OFF to skip ROOT
    cmake --build build -j8
    ctest --test-dir build                              # parity with the Fortran outputs in data/e94010

## Run
    build/radcor unpol --dir data/e94010/unpol --out e94010.dat
    build/radcor unpol --input mydata.csv --precision double --set legacy=0 \
        --set collisional_loss=0 --set finite_acc=0 --details details.txt --root details.root
    build/radcor check unpol --input mydata.csv          # validate an input
    build/radcor convert unpol --input ele-input --to mydata.csv   # ele <-> csv <-> root
    build/radcor -h                                     # all options

Input formats: the Fortran files (ele-input; radcor_input.dat + ds.inp), a CSV table with named
columns, or a ROOT tree; from Python, arrays via python/radcor.py (manual, chapter 3).

## Layout
    apps/radcor.cpp            command-line front end (run, check, convert)
    include/radcor/Input.h     the input, independent of format (RunInput), formats, validation
    include/radcor/Run.h       drivers: RunInput -> engine -> outputs
    include/radcor/Engine.h    the unfolding engine (radiating formula, model, iteration)
    include/radcor/Options.h   every switch and parameter, with defaults
    include/radcor/Interp.h    spline and polynomial interpolation of the output stage
    src/                       formats, drivers, option parsing, ROOT input/output
    python/radcor.py           Python interface
    fortran/                   the original Fortran sources, Makefile
    data/e94010/               E94-010 inputs and the Fortran reference outputs
    examples/                  reading the outputs; the library from C++ and Python
    tests/                     ctest: parity, input formats, Python, Fortran rebuild
    doc/                       the manual; K. Slifer's description of his changes
