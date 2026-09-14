#!/usr/bin/env bash
#
# Tests scripts/tif2lcp.py and the --recompute-lcp-stats safety net.
#
# The strategy is a round-trip: take an .lcp whose FARSITE output is known to
# reproduce the committed references, push it out to GeoTIFF, rebuild an .lcp
# from that with the converter, and require that FARSITE produce the same
# results from the rebuilt file. That checks the raster, the header geometry and
# the statistics block all at once, against a known-good answer.
#
# Skips (exit 0) if GDAL is unavailable.

set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LCP="$REPO/examples/Panther/cust/input/238648_cust_fm94.lcp"
REFDIR="$REPO/examples/Panther/cust/output"
RUNFILE="$REPO/examples/Panther/runPanther.txt"

command -v gdal_translate >/dev/null 2>&1 || { echo "SKIP: gdal_translate not found"; exit 0; }
python3 -c "from osgeo import gdal" >/dev/null 2>&1 || { echo "SKIP: GDAL python bindings not found"; exit 0; }
[ -f "$LCP" ] || { echo "SKIP: $LCP missing"; exit 0; }
[ -x "$REPO/src/TestFARSITE" ] || { echo "SKIP: build TestFARSITE first (make -C src)"; exit 0; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/farsite-tif2lcp.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
fail=0
ok()   { echo "  PASS  $1"; }
bad()  { echo "  FAIL  $1"; fail=1; }

# The example command file hardcodes $scriptRoot; stage a copy pointing at $WORK.
mkdir -p "$WORK/farsite"
cp -R "$REPO/examples" "$WORK/farsite/"
grep -rlF '$scriptRoot' "$WORK/farsite/examples" --exclude-dir=.git \
  | REPL="$WORK" xargs perl -pi -e 's/\$scriptRoot/$ENV{REPL}/g'

# Build a run-command line: same inputs as Panther "cust", but our LCP + outdir.
mk_run() { # $1=lcp  $2=out base  -> prints run file path
  local rf="$WORK/run_$(basename "$2").txt"
  head -1 "$WORK/farsite/examples/Panther/runPanther.txt" \
    | awk -v L="$1" -v O="$2" '{print L" "$2" "$3" "$4" "O" "$6}' > "$rf"
  mkdir -p "$(dirname "$2")"
  echo "$rf"
}

cmp_refs() { # $1=produced prefix  $2=label ; compares all .asc to the references
  local n=0 d=0 f r
  for f in "$1"_*.asc; do
    [ -f "$f" ] || continue
    r="$REFDIR/cust_ignit2_${f##*"$(basename "$1")"_}"
    [ -f "$r" ] || continue
    n=$((n+1)); cmp -s "$f" "$r" || { d=$((d+1)); echo "        differs: $(basename "$r")"; }
  done
  if [ "$n" -eq 0 ]; then bad "$2 (no outputs produced)"
  elif [ "$d" -eq 0 ]; then ok "$2 ($n/$n ASCII grids identical to reference)"
  else bad "$2 ($d of $n ASCII grids differ)"; fi
}

echo "== tif2lcp round-trip =="
gdal_translate -q -of GTiff "$LCP" "$WORK/stack.tif" || { bad "gdal_translate to GeoTIFF"; exit 1; }

if python3 "$REPO/scripts/tif2lcp.py" --stack "$WORK/stack.tif" -o "$WORK/new.lcp" \
     --latitude 37 >"$WORK/conv.log" 2>&1; then
  ok "converter ran"
else
  bad "converter failed:"; sed 's/^/        /' "$WORK/conv.log"; exit 1
fi

# 1. raster payload must survive byte for byte
if [ "$(cmp -l "$WORK/new.lcp" "$LCP" 2>/dev/null | awk '$1>7316' | wc -l)" -eq 0 ]; then
  ok "raster section byte-identical to the original LCP"
else
  bad "raster section differs from the original LCP"
fi

# 2. size must match exactly (header 7316 + ncols*nrows*NumVals*2)
if [ "$(stat -f%z "$WORK/new.lcp" 2>/dev/null || stat -c%s "$WORK/new.lcp")" \
   = "$(stat -f%z "$LCP" 2>/dev/null || stat -c%s "$LCP")" ]; then
  ok "file size matches"
else
  bad "file size differs"
fi

# 3. GDAL must be able to re-read the statistics block we wrote. GDAL's own LCP
#    writer fails this: it returns FUEL_MODEL_VALUES= empty.
vals="$(gdalinfo "$WORK/new.lcp" 2>/dev/null | sed -n 's/ *FUEL_MODEL_VALUES=//p')"
if [ -n "$vals" ] && [ "$vals" = "$(gdalinfo "$LCP" 2>/dev/null | sed -n 's/ *FUEL_MODEL_VALUES=//p')" ]; then
  ok "GDAL reads back the fuel category list, matching the original"
else
  bad "fuel category list not readable/mismatched (got '${vals:-<empty>}')"
fi

# 4. the real test: FARSITE output from the rebuilt LCP
rf="$(mk_run "$WORK/new.lcp" "$WORK/out_conv/cv")"
(cd "$REPO/src" && ./TestFARSITE "$rf") >"$WORK/run_conv.log" 2>&1 \
  || bad "FARSITE run on converted LCP failed"
cmp_refs "$WORK/out_conv/cv" "FARSITE output from converted LCP"

echo
echo "== --recompute-lcp-stats safety net =="
# Corrupt only the statistics block (bytes 44..4164), leaving the raster intact.
# This is what GDAL's LCP writer effectively does.
python3 - "$LCP" "$WORK/corrupt.lcp" <<'PY'
import sys
data = bytearray(open(sys.argv[1], 'rb').read())
for i in range(44, 4164):
    data[i] = 0
open(sys.argv[2], 'wb').write(bytes(data))
PY

rf="$(mk_run "$WORK/corrupt.lcp" "$WORK/out_bad/bd")"
(cd "$REPO/src" && ./TestFARSITE "$rf") >"$WORK/run_bad.log" 2>&1
bad_differs=0
for f in "$WORK"/out_bad/bd_*.asc; do
  [ -f "$f" ] || continue
  r="$REFDIR/cust_ignit2_${f##*bd_}"
  [ -f "$r" ] || continue
  cmp -s "$f" "$r" || bad_differs=$((bad_differs+1))
done
if [ "$bad_differs" -gt 0 ]; then
  ok "corrupt statistics block does change results ($bad_differs grids) -- the test is meaningful"
else
  bad "corrupt statistics block changed nothing; this test proves nothing"
fi

rf="$(mk_run "$WORK/corrupt.lcp" "$WORK/out_fix/fx")"
(cd "$REPO/src" && ./TestFARSITE "$rf" --recompute-lcp-stats) >"$WORK/run_fix.log" 2>&1 \
  || bad "FARSITE run with --recompute-lcp-stats failed"
grep -q "Recomputing landscape statistics" "$WORK/run_fix.log" \
  && ok "flag was honoured" || bad "flag had no effect (no log line)"
cmp_refs "$WORK/out_fix/fx" "recovered output with --recompute-lcp-stats"

echo
echo "== oracle: the flag must be a no-op on a well-formed LCP =="
rf="$(mk_run "$LCP" "$WORK/out_oracle/oc")"
(cd "$REPO/src" && ./TestFARSITE "$rf" --recompute-lcp-stats) >"$WORK/run_oracle.log" 2>&1 \
  || bad "oracle run failed"
cmp_refs "$WORK/out_oracle/oc" "recomputed stats agree with a good header"

echo
[ $fail -eq 0 ] && echo "TIF2LCP TESTS PASS" || echo "TIF2LCP TESTS FAIL"
exit $fail
