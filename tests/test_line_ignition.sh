#!/usr/bin/env bash
#
# Guards the polyline (line-source) ignition path in IgnitionFile::ShapeInput.
#
# That branch builds a degenerate "out and back" perimeter: it writes the N
# vertices forward into slots 0..N-1, then writes the N-2 interior vertices
# backward into slots N..2N-3, for 2N-2 points total. The backward write used
# to be guarded by `k < end - 1`, but `k` is only ever assigned inside the
# Point branch, so for a polyline-only shapefile it was read uninitialised.
# Half the perimeter was then left unwritten while SetNumPoints claimed all of
# it, and the simulation never terminated. Fixed to `j < end - 1`, matching the
# equivalent guard in ArcLine().
#
# The test is differential: a 2-vertex line goes through a completely separate
# code path that never had the bug, so describing the same physical line with 2
# and with 3 collinear vertices must produce the same fire.
#
# Skips (exit 0) if GDAL or the binary is missing.

set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LCP="$REPO/examples/Panther/cust/input/238648_cust_fm94.lcp"
INPUT="$REPO/examples/Panther/cust/input/cust.input"
IGN="$REPO/examples/Panther/cust/input/ignit2.shp"

python3 -c "from osgeo import ogr" >/dev/null 2>&1 || { echo "SKIP: GDAL python bindings not found"; exit 0; }
[ -x "$REPO/src/TestFARSITE" ] || { echo "SKIP: build TestFARSITE first (make -C src)"; exit 0; }
for f in "$LCP" "$INPUT" "$IGN"; do
  [ -f "$f" ] || { echo "SKIP: missing $f"; exit 0; }
done

WORK="$(mktemp -d "${TMPDIR:-/tmp}/farsite-lineign.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
fail=0

# The example .input references its RAWS/wind files by relative or $scriptRoot
# path, so stage the tree the same way the other tests do.
mkdir -p "$WORK/farsite"
cp -R "$REPO/examples" "$WORK/farsite/"
grep -rlF '$scriptRoot' "$WORK/farsite/examples" --exclude-dir=.git 2>/dev/null \
  | REPL="$WORK" xargs -r perl -pi -e 's/\$scriptRoot/$ENV{REPL}/g' 2>/dev/null || true
W_INPUT="$WORK/farsite/examples/Panther/cust/input/cust.input"
W_LCP="$WORK/farsite/examples/Panther/cust/input/238648_cust_fm94.lcp"

# Build collinear lines of 2 and 3 vertices describing one 600 m line, centred
# on the existing example ignition so it lands on burnable fuel.
python3 - "$WORK" "$WORK/farsite/examples/Panther/cust/input/ignit2.shp" <<'PY'
import contextlib, io, sys
from osgeo import ogr, osr
with contextlib.redirect_stderr(io.StringIO()):
    ogr.UseExceptions()
work, ign = sys.argv[1], sys.argv[2]
ds = ogr.Open(ign); lyr = ds.GetLayer(0)
xs, ys, n = 0.0, 0.0, 0
for feat in lyr:
    g = feat.GetGeometryRef()
    c = g.Centroid()
    xs += c.GetX(); ys += c.GetY(); n += 1
cx, cy = xs / n, ys / n
srs = lyr.GetSpatialRef()
ds = None
for nv in (2, 3):
    out = "%s/line%d.shp" % (work, nv)
    d = ogr.GetDriverByName("ESRI Shapefile").CreateDataSource(out)
    l = d.CreateLayer("l", srs, ogr.wkbLineString)
    l.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))
    g = ogr.Geometry(ogr.wkbLineString)
    for i in range(nv):
        g.AddPoint_2D(cx - 300.0 + 600.0 * i / (nv - 1), cy)
    f = ogr.Feature(l.GetLayerDefn()); f.SetGeometry(g); f.SetField("id", 0)
    l.CreateFeature(f); d = None
print("  ignition lines centred on (%.1f, %.1f)" % (cx, cy))
PY

echo "== line-source ignition =="
for nv in 2 3; do
  mkdir -p "$WORK/out$nv"
  echo "$W_LCP $W_INPUT $WORK/line$nv.shp 0 $WORK/out$nv/l 1" > "$WORK/run$nv.txt"
  # A regression here manifests as a run that never ends, so bound it.
  ( cd "$REPO/src" && ./TestFARSITE "$WORK/run$nv.txt" ) >"$WORK/run$nv.log" 2>&1 &
  pid=$!
  waited=0
  while kill -0 "$pid" 2>/dev/null; do
    sleep 2; waited=$((waited+2))
    if [ $waited -ge 120 ]; then
      kill -9 "$pid" 2>/dev/null
      echo "  FAIL  $nv-vertex line did not finish in 120s (the uninitialised-\`k\` symptom)"
      fail=1; break
    fi
  done
  wait "$pid" 2>/dev/null
  if [ $waited -lt 120 ]; then
    if [ -f "$WORK/out$nv/l_ArrivalTime.asc" ]; then
      echo "  PASS  $nv-vertex line completed in ~${waited}s and wrote output"
    else
      echo "  FAIL  $nv-vertex line produced no ArrivalTime grid"; fail=1
    fi
  fi
done

if [ -f "$WORK/out2/l_ArrivalTime.asc" ] && [ -f "$WORK/out3/l_ArrivalTime.asc" ]; then
  python3 - "$WORK/out2/l_ArrivalTime.asc" "$WORK/out3/l_ArrivalTime.asc" <<'PY'
import sys
def cells(path):
    t = open(path).read().split()
    h = {t[i].lower(): t[i+1] for i in range(0, 12, 2)}
    nc = int(h['ncols']); nod = float(h['nodata_value'])
    v = [float(x) for x in t[12:]]
    return {i for i, x in enumerate(v) if x != nod and x > 0}
a, b = cells(sys.argv[1]), cells(sys.argv[2])
if not a or not b:
    print("  FAIL  one of the runs burned nothing (%d vs %d cells)" % (len(a), len(b)))
    sys.exit(1)
j = len(a & b) / float(len(a | b))
print("  %s  2-vertex vs 3-vertex agreement: Jaccard %.4f (%d vs %d cells)"
      % ("PASS" if j > 0.99 else "FAIL", j, len(a), len(b)))
sys.exit(0 if j > 0.99 else 1)
PY
  [ $? -eq 0 ] || fail=1
fi

echo
[ $fail -eq 0 ] && echo "LINE IGNITION TESTS PASS" || echo "LINE IGNITION TESTS FAIL"
exit $fail
