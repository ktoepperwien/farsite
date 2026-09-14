#!/usr/bin/env bash
#
# Regression test for FARSITE.
#
# The repo ships reference outputs next to each example, generated on x86 Linux
# with gcc at -O0. Re-running them elsewhere splits into two groups:
#
#   STRICT cases   reproduce the references bit-for-bit. They are the real
#                  regression gate: any difference in a .asc grid is a bug.
#
#   SENSITIVE cases do not, and cannot, reproduce bit-for-bit off the original
#                  platform. FARSITE's spread front is threshold-driven, so
#                  last-bit floating-point differences get amplified into
#                  macroscopic ones. Measured on Apple silicon: Panther's
#                  test1177973 differs from the references in 9 of 10 grids,
#                  and it differs that much between -O0 and -O2 *on the same
#                  machine with identical sources*. These cases are therefore
#                  run as smoke tests: they must complete and produce output,
#                  and their divergence is reported for eyeballing, but they do
#                  not fail the run. Set STRICT_ALL=1 to gate on them anyway
#                  (only meaningful on x86 Linux/gcc -O0).
#
# .fbg files are raw float32 and .shp files store doubles, so both are compared
# numerically rather than byte-wise: even STRICT cases legitimately differ in
# the last bit on a few cells. .asc/.csv/.dbf are text and are compared exactly.
# .shx is skipped: it carries no data beyond the .shp bbox and record offsets.
#
# Usage:  tests/run_regression.sh [case ...]

set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/farsite-regress.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

FBG_RTOL="${FBG_RTOL:-1e-5}"
SHP_RTOL="${SHP_RTOL:-1e-6}"
STRICT_ALL="${STRICT_ALL:-0}"

# name : strictness : run-command file (relative to examples/)
CASES=(
  "panther:strict:Panther/runPanther.txt"
  "cougarCreek:sensitive:cougarCreek/run_cougarCreek_polygonIgnit.txt"
  "flatland:sensitive:flatland/caseRunFiles/centerIgnit/b_RAWS/singleWinds/run_flatland_centerIgnit_3mph0deg7hr_RAWS.txt"
)

# Output basenames whose references are known not to reproduce off x86/gcc -O0.
# Everything else in a "strict" case must match exactly.
SENSITIVE_PREFIXES="test1177973_NTFBIgnition"

fail=0; n_exact=0; n_tol=0; n_soft=0; n_hard=0

echo "== building =="
make -C "$REPO/src" -j"$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)" >"$WORK/build.log" 2>&1 || {
  echo "BUILD FAILED:"; tail -30 "$WORK/build.log"; exit 1; }
echo "ok: $REPO/src/TestFARSITE"

echo "== staging examples into $WORK =="
mkdir -p "$WORK/farsite"
cp -R "$REPO/examples" "$WORK/farsite/"
# grep -F is required: without it, BSD grep reads the '$' as an end-of-line anchor.
grep -rlF '$scriptRoot' "$WORK/farsite/examples" --exclude-dir=.git \
  | REPL="$WORK" xargs perl -pi -e 's/\$scriptRoot/$ENV{REPL}/g'
STAMP="$WORK/.stamp"; touch "$STAMP"

compare_fbg() {
  python3 - "$1" "$2" "$FBG_RTOL" <<'PY'
import struct, sys
new, ref, rtol = sys.argv[1], sys.argv[2], float(sys.argv[3])
a, b = open(new,'rb').read(), open(ref,'rb').read()
if len(a) != len(b):
    print("size differs: %d vs %d" % (len(a), len(b))); sys.exit(1)
off = 32
n = (len(a)-off)//4
A = struct.unpack('<%df'%n, a[off:off+4*n]); B = struct.unpack('<%df'%n, b[off:off+4*n])
worst = ndiff = 0
for x,y in zip(A,B):
    if abs(x) > 1e30 or abs(y) > 1e30: continue
    if x != y:
        ndiff += 1; worst = max(worst, abs(x-y)/max(abs(y),1e-9))
if worst > rtol:
    print("%d cells differ, max rel %.3g > tol %.3g" % (ndiff, worst, rtol)); sys.exit(1)
if ndiff: print("%d/%d cells differ, max rel %.3g (within tol)" % (ndiff, n, worst))
sys.exit(0)
PY
}

is_sensitive_file() {
  for p in $SENSITIVE_PREFIXES; do
    case "$(basename "$1")" in $p*) return 0;; esac
  done
  return 1
}

for entry in "${CASES[@]}"; do
  name="${entry%%:*}"; rest="${entry#*:}"
  strictness="${rest%%:*}"; runfile="${rest#*:}"
  if [ $# -gt 0 ]; then
    want=0; for a in "$@"; do [ "$a" = "$name" ] && want=1; done
    [ $want -eq 1 ] || continue
  fi

  echo; echo "== case: $name ($strictness) =="
  cmd="$WORK/farsite/examples/$runfile"
  [ -f "$cmd" ] || { echo "SKIP: no run file at $runfile"; continue; }
  while read -r base; do [ -n "$base" ] && mkdir -p "$(dirname "$base")"; done < <(awk '{print $5}' "$cmd")

  if ! (cd "$REPO/src" && ./TestFARSITE "$cmd") >"$WORK/$name.log" 2>&1; then
    echo "  RUN FAILED (see $WORK/$name.log)"; tail -20 "$WORK/$name.log"; fail=1; n_hard=$((n_hard+1)); continue
  fi

  produced_any=0
  while read -r produced; do
    produced_any=1
    rel="${produced#$WORK/farsite/}"; ref="$REPO/$rel"
    [ -f "$ref" ] || continue

    # decide whether a mismatch here should fail the run
    soft=0
    [ "$strictness" = "sensitive" ] && soft=1
    is_sensitive_file "$produced" && soft=1
    [ "$STRICT_ALL" = "1" ] && soft=0

    case "$produced" in
      *.fbg)
        if out=$(compare_fbg "$produced" "$ref"); then
          if [ -n "$out" ]; then n_tol=$((n_tol+1)); echo "  ~ $(basename "$rel"): $out"; else n_exact=$((n_exact+1)); fi
        elif [ $soft -eq 1 ]; then n_soft=$((n_soft+1)); echo "  DIVERGED (not gating) $(basename "$rel"): $out"
        else n_hard=$((n_hard+1)); fail=1; echo "  FAIL $(basename "$rel"): $out"
        fi ;;
      *.shp)
        if out=$(python3 "$REPO/tests/compare_shp.py" "$produced" "$ref" "$SHP_RTOL" 2>&1); then
          if [ -n "$out" ]; then n_tol=$((n_tol+1)); echo "  ~ $(basename "$rel"): $out"; else n_exact=$((n_exact+1)); fi
        elif [ $soft -eq 1 ]; then n_soft=$((n_soft+1)); echo "  DIVERGED (not gating) $(basename "$rel"): $out"
        else n_hard=$((n_hard+1)); fail=1; echo "  FAIL $(basename "$rel"): $out"
        fi ;;
      *.shx)
        ;;
      *)
        if cmp -s "$produced" "$ref"; then n_exact=$((n_exact+1))
        elif [ $soft -eq 1 ]; then n_soft=$((n_soft+1)); echo "  DIVERGED (not gating) $(basename "$rel")"
        else n_hard=$((n_hard+1)); fail=1; echo "  FAIL $(basename "$rel") (not byte-identical)"
        fi ;;
    esac
  done < <(find "$WORK/farsite/examples" -newer "$STAMP" -type f \
             \( -name '*.asc' -o -name '*.fbg' -o -name '*.csv' -o -name '*.shp' -o -name '*.dbf' -o -name '*.shx' \) | sort)
  [ $produced_any -eq 1 ] || { echo "  FAIL: run produced no output files"; fail=1; n_hard=$((n_hard+1)); }
done

echo
echo "==================================================="
printf "bit-identical: %d   within-tolerance: %d\n" "$n_exact" "$n_tol"
printf "diverged (platform-sensitive, not gating): %d\n" "$n_soft"
printf "failures: %d\n" "$n_hard"
[ $fail -eq 0 ] && echo "REGRESSION PASS" || echo "REGRESSION FAIL"
exit $fail
