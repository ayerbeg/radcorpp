#!/usr/bin/env bash
# Parity test: the C++ in float mode with the default (legacy) settings must reproduce the
# Fortran output byte for byte on the E94-010 inputs shipped in data/e94010.
# Usage: tests/run_parity.sh   (from anywhere; needs build/radcor)
set -euo pipefail
here=$(cd "$(dirname "$0")/.." && pwd)
data=$here/data/e94010
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
"$here/build/radcor" unpol --dir "$data/unpol" --out "$tmp/u.dat" -q
"$here/build/radcor" pol   --dir "$data/pol"   --out "$tmp/p.dat" -q
fail=0
cmp -s "$tmp/u.dat" "$data/unpol/expected/E94010_xs_unfolded.dat" \
  && echo "PASS unpolarized: byte-identical to Fortran" || { echo "FAIL unpolarized"; fail=1; }
cmp -s "$tmp/p.dat" "$data/pol/expected/ds.rc3" \
  && echo "PASS polarized: byte-identical to Fortran" || { echo "FAIL polarized"; fail=1; }
exit $fail
