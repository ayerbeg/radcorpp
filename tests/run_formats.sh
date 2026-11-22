#!/usr/bin/env bash
# Input formats: convert the E94-010 inputs ele/pol -> csv -> root -> csv and back to the Fortran
# formats; every conversion must be lossless and every format must give the Fortran output.
set -euo pipefail
here=$(cd "$(dirname "$0")/.." && pwd); R=$here/build/radcor; data=$here/data/e94010
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
fail=0
ok() { echo "PASS $1"; }; ko() { echo "FAIL $1"; fail=1; }
# unpolarized
"$R" convert unpol --dir "$data/unpol" --to "$tmp/u.csv" -q
"$R" convert unpol --input "$tmp/u.csv" --to "$tmp/u.root" -q 2>/dev/null || true
"$R" convert unpol --input "$tmp/u.csv" --to "$tmp/ele" -q
"$R" convert unpol --input "$tmp/ele" --to "$tmp/u2.csv" -q
cmp -s "$tmp/u.csv" "$tmp/u2.csv" && ok "unpol ele->csv->ele->csv lossless" || ko "unpol ele->csv->ele->csv"
mkdir -p "$tmp/u/output"; cp -r "$data/unpol/input" "$tmp/u/"; cp "$tmp/u.csv" "$tmp/u/"
"$R" unpol --dir "$tmp/u" --input u.csv --out a.dat -q
cmp -s "$tmp/u/a.dat" "$data/unpol/expected/E94010_xs_unfolded.dat" && ok "unpol csv input == Fortran" || ko "unpol csv input"
if [ -f "$tmp/u.root" ]; then
  "$R" convert unpol --input "$tmp/u.root" --to "$tmp/u3.csv" -q
  cmp -s "$tmp/u.csv" "$tmp/u3.csv" && ok "unpol csv->root->csv lossless" || ko "unpol csv->root->csv"
  cp "$tmp/u.root" "$tmp/u/"; "$R" unpol --dir "$tmp/u" --input u.root --out b.dat -q
  cmp -s "$tmp/u/b.dat" "$data/unpol/expected/E94010_xs_unfolded.dat" && ok "unpol root input == Fortran" || ko "unpol root input"
else echo "SKIP root (built without ROOT)"; fi
# polarized
"$R" convert pol --dir "$data/pol" --to "$tmp/p.csv" -q
"$R" convert pol --input "$tmp/p.csv" --to "$tmp/hdr" --to-data "$tmp/ds" -q
"$R" convert pol --input "$tmp/hdr" --data "$tmp/ds" --to "$tmp/p2.csv" -q
cmp -s "$tmp/p.csv" "$tmp/p2.csv" && ok "pol pol->csv->pol->csv lossless" || ko "pol pol->csv->pol->csv"
"$R" pol --input "$tmp/p.csv" --out "$tmp/p.dat" -q
cmp -s "$tmp/p.dat" "$data/pol/expected/ds.rc3" && ok "pol csv input == Fortran" || ko "pol csv input"
# validation must reject a reordered input
{ sed -n 1,3p "$tmp/u.csv"; sed -n '4,13p' "$tmp/u.csv" | sort -t, -k3 -r; sed -n '14,$p' "$tmp/u.csv"; } > "$tmp/bad.csv"
if "$R" check unpol --input "$tmp/bad.csv" > /dev/null; then ko "validation accepted unsorted nu"; else ok "validation rejects unsorted nu"; fi
exit $fail
