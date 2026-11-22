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
    build/radcor unpol --dir mydata --precision double --set legacy=0 \
        --set collisional_loss=0 --set finite_acc=0 --details details.txt --root details.root
    build/radcor -h                                     # all options

## Layout
    apps/radcor.cpp            command-line front end
    include/radcor/Engine.h    the unfolding engine (radiating formula, model, iteration)
    include/radcor/Options.h   every switch and parameter, with defaults
    include/radcor/Interp.h    spline and polynomial interpolation of the output stage
    include/radcor/Legacy.h    drivers that read and write the Fortran file formats
    src/                       drivers, option parsing, details and ROOT writers
    data/e94010/               E94-010 inputs and the Fortran reference outputs
    examples/                  reading the outputs; using the engine from your own program
    tests/                     parity test and comparison scripts
    doc/radcor_manual.pdf      the manual
