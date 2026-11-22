#!/usr/bin/env bash
# Build the Fortran programs from fortran/, run them on the E94-010 inputs, and check that
# (1) they reproduce data/e94010/*/expected and (2) the C++ float mode reproduces them.
# Needs gfortran and build/radcor.
set -euo pipefail
here=$(cd "$(dirname "$0")/.." && pwd)
data=$here/data/e94010
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
make -s -C "$here/fortran" FC="${FC:-gfortran}" >/dev/null 2>&1 || make -C "$here/fortran"
mkdir -p "$tmp/u/output" "$tmp/p/output"
cp -r "$data/unpol/ele-input" "$data/unpol/input" "$tmp/u/"
cp -r "$data/pol/input" "$tmp/p/"
( cd "$tmp/u" && "$here/fortran/unpolarized/nradcor" > /dev/null )
( cd "$tmp/p" && "$here/fortran/polarized/nradcor" > /dev/null )
"$here/build/radcor" unpol --dir "$tmp/u" --out cpp.dat -q
"$here/build/radcor" pol   --dir "$tmp/p" --out cpp.dat -q
fail=0
check() { cmp -s "$1" "$2" && echo "PASS $3" || { echo "FAIL $3"; fail=1; }; }
check "$tmp/u/output/E94010_xs_unfolded.dat" "$data/unpol/expected/E94010_xs_unfolded.dat" "Fortran U == expected"
check "$tmp/p/output/ds.rc3" "$data/pol/expected/ds.rc3" "Fortran P == expected"
check "$tmp/u/cpp.dat" "$tmp/u/output/E94010_xs_unfolded.dat" "C++ U == Fortran U"
check "$tmp/p/cpp.dat" "$tmp/p/output/ds.rc3" "C++ P == Fortran P"
exit $fail
