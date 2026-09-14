#!/usr/bin/env python3
"""Create or validate a FARSITE ignition shapefile against a landscape file.

Two modes:

  circle    (default) write a polygon circle, by default at the centre of the
            landscape with a 200 m radius
  check     read an existing shapefile, verify it matches what FARSITE's reader
            accepts, and write a cleaned copy

Both modes report where the ignition sits relative to the landscape edge and
what fuel lies underneath it, because the two ways an ignition silently fails
are landing outside the landscape (everything off-grid reads as fuel -9999,
which converts to unburnable) and landing entirely on a non-burnable fuel model.

Examples
--------
  # 200 m circle at the centre of the landscape
  python3 scripts/make_ignition.py --lcp examples/Dish/0_Dish/0_Dish.lcp \
      -o examples/Dish/0_Dish/ignition.shp

  # 500 m circle at a chosen point
  python3 scripts/make_ignition.py --lcp case.lcp -o ign.shp \
      --center -2266455 1910685 --radius 500

  # validate and clean something you already have
  python3 scripts/make_ignition.py --lcp case.lcp -o ign_clean.shp \
      --mode check --input my_perimeter.shp

What FARSITE accepts
--------------------
`IgnitionFile::ShapeInput` switches on the shapefile's type:

  Polygon (incl. Z/M)  used directly as an area ignition. Winding is normalised
                       internally: FARSITE computes a signed area and reverses
                       the ring when it comes out negative, so it wants
                       counter-clockwise but accepts either.
  PolyLine / Arc       treated as a line source: the fire spreads outward
                       from the line. Requires the fsxwignt.cpp:191 fix.
  Point / MultiPoint   each point becomes a 10-vertex circle of radius
                       `startsize`, which is hard-coded to 1 (one map unit).
                       On a 30 m landscape that is far smaller than a cell, so
                       point ignitions are near-degenerate -- prefer a polygon.
                       Use --point-to-circle to convert them here instead.

FARSITE removes duplicate vertices and redistributes vertex spacing itself
(RemoveIdenticalPoints, DensityControl), so neither needs to be exact here.
"""

import argparse
import math
import os
import struct
import sys

try:
    import contextlib
    import io
    from osgeo import ogr, osr
    with contextlib.redirect_stderr(io.StringIO()):
        ogr.UseExceptions()
        osr.UseExceptions()
except ImportError:
    sys.exit("error: GDAL/OGR Python bindings are required (pip install gdal, or conda install gdal)")

LCP_HEADER_SIZE = 7316
LCP_OFF_CROWN = 0
LCP_OFF_NUMEAST = 4164
LCP_OFF_EASTUTM = 4172

# Fuel model codes FARSITE cannot carry fire. 91-99 are the standard
# non-burnable block (urban, snow, agriculture, water, barren); 0 and -9999
# are no-data. GetFuelConversion maps anything outside [0,257) to -1.
NONBURNABLE = set(range(91, 100)) | {0, -9999, -1}


class Lcp(object):
    """Just enough of the landscape file to place and check an ignition."""

    def __init__(self, path):
        self.path = path
        with open(path, "rb") as fh:
            head = fh.read(LCP_HEADER_SIZE)
        if len(head) < LCP_HEADER_SIZE:
            sys.exit("error: %s is too short to be an .lcp" % path)
        crown, ground = struct.unpack_from("<2i", head, LCP_OFF_CROWN)
        self.numeast, self.numnorth = struct.unpack_from("<2i", head, LCP_OFF_NUMEAST)
        self.east, self.west, self.north, self.south = struct.unpack_from(
            "<4d", head, LCP_OFF_EASTUTM)
        if self.numeast <= 0 or self.numnorth <= 0:
            sys.exit("error: %s reports a %dx%d grid; not a readable .lcp"
                     % (path, self.numeast, self.numnorth))
        # Matches ReadHeader: resolution comes from the corners, not XResol.
        self.resx = (self.east - self.west) / self.numeast
        self.resy = (self.north - self.south) / self.numnorth
        self.numvals = {(20, 20): 5, (20, 21): 7, (21, 20): 8, (21, 21): 10}.get(
            (crown, ground))
        if self.numvals is None:
            sys.exit("error: %s has CrownFuels=%d GroundFuels=%d, expected 20/21"
                     % (path, crown, ground))
        self._fuel_cache = {}

    @property
    def center(self):
        return (self.west + self.east) / 2.0, (self.south + self.north) / 2.0

    def contains(self, x, y):
        return self.west <= x <= self.east and self.south <= y <= self.north

    def edge_margin(self, x, y):
        """Distance from (x, y) to the nearest landscape edge; negative if outside."""
        return min(x - self.west, self.east - x, y - self.south, self.north - y)

    def cell_of(self, x, y):
        col = int((x - self.west) / self.resx)
        row = int((self.north - y) / self.resy)
        return col, row

    def fuel_at_cell(self, col, row):
        """Fuel model at a cell, read straight from the interleaved raster."""
        if not (0 <= col < self.numeast and 0 <= row < self.numnorth):
            return None
        if row not in self._fuel_cache:
            with open(self.path, "rb") as fh:
                fh.seek(LCP_HEADER_SIZE + row * self.numeast * self.numvals * 2)
                raw = fh.read(self.numeast * self.numvals * 2)
            vals = struct.unpack("<%dh" % (self.numeast * self.numvals), raw)
            self._fuel_cache[row] = vals[3::self.numvals]   # index 3 = fuel model
        return self._fuel_cache[row][col]

    def describe(self):
        print("landscape: %s" % self.path)
        print("  %d x %d cells, res %g x %g" % (self.numeast, self.numnorth, self.resx, self.resy))
        print("  extent x[%.1f, %.1f] y[%.1f, %.1f]" % (self.west, self.east, self.south, self.north))
        print("  centre (%.1f, %.1f)" % self.center)


def find_crs(lcp_path, explicit=None):
    """Best-effort CRS for the output .prj.

    The .lcp format stores no CRS at all, so look for one alongside it: the
    crop_info.json that crop_landfire.py writes, or a sibling .prj.
    """
    if explicit:
        srs = osr.SpatialReference()
        if srs.SetFromUserInput(explicit) != 0:
            sys.exit("error: could not interpret --crs %r" % explicit)
        return srs, "--crs"
    folder = os.path.dirname(os.path.abspath(lcp_path))
    cand = os.path.join(folder, "crop_info.json")
    if os.path.exists(cand):
        try:
            import json
            wkt = json.load(open(cand)).get("crs_wkt")
            if wkt:
                return osr.SpatialReference(wkt=wkt), os.path.basename(cand)
        except Exception:
            pass
    cand = os.path.splitext(lcp_path)[0] + ".prj"
    if os.path.exists(cand):
        srs = osr.SpatialReference()
        if srs.SetFromUserInput(open(cand).read()) == 0:
            return srs, os.path.basename(cand)
    return None, None


def signed_area(pts):
    """Shoelace area. Positive is counter-clockwise, which is what arp() wants."""
    a = 0.0
    for i in range(len(pts) - 1):
        a += pts[i][0] * pts[i + 1][1] - pts[i + 1][0] * pts[i][1]
    return a / 2.0


def perimeter(pts):
    return sum(math.hypot(pts[i + 1][0] - pts[i][0], pts[i + 1][1] - pts[i][1])
               for i in range(len(pts) - 1))


def close_ring(pts):
    return pts if pts and pts[0] == pts[-1] else pts + [pts[0]]


def dedupe(pts, tol=1e-9):
    out = [pts[0]]
    for p in pts[1:]:
        if math.hypot(p[0] - out[-1][0], p[1] - out[-1][1]) > tol:
            out.append(p)
    return out


def make_circle(cx, cy, radius, nvert, ccw=True):
    step = 2.0 * math.pi / nvert
    pts = [(cx + radius * math.cos(i * step), cy + radius * math.sin(i * step))
           for i in range(nvert)]
    if not ccw:
        pts.reverse()
    return close_ring(pts)


def is_closed(pts, tol=1e-9):
    """True if the vertex list forms a closed ring (a polygon, not a line)."""
    return (len(pts) >= 4
            and math.hypot(pts[0][0] - pts[-1][0], pts[0][1] - pts[-1][1]) <= tol)


def cells_covered(lcp, geoms):
    """Set of (col, row) the ignition occupies.

    Closed rings use a point-in-polygon test on cell centres. Open geometry --
    a line source, or a bare point -- has no area, so instead walk its segments
    and take every cell they pass through; testing cell centres against a line
    would report nothing covered.
    """
    cells = set()
    step = min(lcp.resx, lcp.resy) / 3.0
    for pts in geoms:
        if is_closed(pts):
            xs = [p[0] for p in pts]
            ys = [p[1] for p in pts]
            c0, r0 = lcp.cell_of(min(xs), max(ys))
            c1, r1 = lcp.cell_of(max(xs), min(ys))
            for row in range(max(0, r0), min(lcp.numnorth, r1 + 1)):
                for col in range(max(0, c0), min(lcp.numeast, c1 + 1)):
                    x = lcp.west + (col + 0.5) * lcp.resx
                    y = lcp.north - (row + 0.5) * lcp.resy
                    if point_in_ring(x, y, pts):
                        cells.add((col, row))
        elif len(pts) == 1:
            col, row = lcp.cell_of(*pts[0])
            if 0 <= col < lcp.numeast and 0 <= row < lcp.numnorth:
                cells.add((col, row))
        else:
            for i in range(len(pts) - 1):
                x0, y0 = pts[i]
                x1, y1 = pts[i + 1]
                seg = math.hypot(x1 - x0, y1 - y0)
                n = max(1, int(math.ceil(seg / step)))
                for k in range(n + 1):
                    t = k / float(n)
                    col, row = lcp.cell_of(x0 + (x1 - x0) * t, y0 + (y1 - y0) * t)
                    if 0 <= col < lcp.numeast and 0 <= row < lcp.numnorth:
                        cells.add((col, row))
    return cells


def report_footprint(lcp, rings, label="ignition"):
    """Edge clearance and underlying fuel for a set of rings. Returns ok/warnings."""
    xs = [p[0] for r in rings for p in r]
    ys = [p[1] for r in rings for p in r]
    minx, maxx, miny, maxy = min(xs), max(xs), min(ys), max(ys)
    print("\n%s footprint:" % label)
    print("  bbox x[%.1f, %.1f] y[%.1f, %.1f]" % (minx, maxx, miny, maxy))

    problems, warnings = [], []
    corners = [(minx, miny), (minx, maxy), (maxx, miny), (maxx, maxy)]
    margins = [lcp.edge_margin(x, y) for x, y in corners]
    margin = min(margins)
    if margin < 0:
        problems.append("extends outside the landscape by %.1f m; everything off-grid "
                        "reads as fuel -9999 and will not burn" % -margin)
    print("  clearance to nearest landscape edge: %.1f m (%.1f cells)"
          % (margin, margin / lcp.resx))

    # Fuel under the footprint, sampled on the landscape grid.
    counts, burnable = {}, 0
    cells = cells_covered(lcp, rings)
    inside_cells = len(cells)
    for col, row in cells:
        f = lcp.fuel_at_cell(col, row)
        counts[f] = counts.get(f, 0) + 1
        if f not in NONBURNABLE:
            burnable += 1
    if inside_cells == 0:
        kind = "area" if any(is_closed(r) for r in rings) else "line/point"
        problems.append("covers no landscape cell (%s ignition, cell is %g m) -- "
                        "too small or entirely off-grid" % (kind, lcp.resx))
    else:
        print("  covers %d landscape cells, %d of them burnable" % (inside_cells, burnable))
        top = sorted(counts.items(), key=lambda kv: -kv[1])[:8]
        print("  fuel models: %s" % ", ".join(
            "%s%s x%d" % (f, "*" if f in NONBURNABLE else "", n) for f, n in top))
        if burnable == 0:
            problems.append("every cell under the ignition is a non-burnable fuel "
                            "model (marked * above); the fire cannot start")
        elif burnable < inside_cells:
            warnings.append("%d of %d cells under the ignition are non-burnable"
                            % (inside_cells - burnable, inside_cells))
    return problems, warnings


def point_in_ring(x, y, ring):
    """Standard ray-casting test; ring is a closed list of (x, y)."""
    inside = False
    n = len(ring) - 1
    for i in range(n):
        x0, y0 = ring[i]
        x1, y1 = ring[i + 1]
        if (y0 > y) != (y1 > y):
            xint = x0 + (y - y0) / (y1 - y0) * (x1 - x0)
            if x < xint:
                inside = not inside
    return inside


def write_polygons(path, rings, srs):
    folder = os.path.dirname(os.path.abspath(path))
    if folder and not os.path.isdir(folder):
        os.makedirs(folder)
    drv = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(path):
        drv.DeleteDataSource(path)
    ds = drv.CreateDataSource(path)
    lyr = ds.CreateLayer(os.path.splitext(os.path.basename(path))[0], srs, ogr.wkbPolygon)
    lyr.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))
    for i, ring_pts in enumerate(rings):
        ring = ogr.Geometry(ogr.wkbLinearRing)
        for x, y in ring_pts:
            ring.AddPoint_2D(x, y)
        poly = ogr.Geometry(ogr.wkbPolygon)
        poly.AddGeometry(ring)
        feat = ogr.Feature(lyr.GetLayerDefn())
        feat.SetGeometry(poly)
        feat.SetField("id", i)
        lyr.CreateFeature(feat)
        feat = None
    ds = None


def write_lines(path, lines, srs):
    folder = os.path.dirname(os.path.abspath(path))
    if folder and not os.path.isdir(folder):
        os.makedirs(folder)
    drv = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(path):
        drv.DeleteDataSource(path)
    ds = drv.CreateDataSource(path)
    lyr = ds.CreateLayer(os.path.splitext(os.path.basename(path))[0], srs,
                         ogr.wkbLineString)
    lyr.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))
    for i, pts in enumerate(lines):
        geom = ogr.Geometry(ogr.wkbLineString)
        for x, y in pts:
            geom.AddPoint_2D(x, y)
        feat = ogr.Feature(lyr.GetLayerDefn())
        feat.SetGeometry(geom)
        feat.SetField("id", i)
        lyr.CreateFeature(feat)
        feat = None
    ds = None


ACCEPTED = {
    ogr.wkbPoint: "Point", ogr.wkbMultiPoint: "MultiPoint",
    ogr.wkbLineString: "PolyLine", ogr.wkbMultiLineString: "PolyLine",
    ogr.wkbPolygon: "Polygon", ogr.wkbMultiPolygon: "Polygon",
}


def read_and_clean(path, lcp, args):
    """Pull rings out of an existing shapefile, normalising to FARSITE's liking."""
    ds = ogr.Open(path)
    if ds is None:
        sys.exit("error: cannot open %s" % path)
    lyr = ds.GetLayer(0)
    gtype = lyr.GetGeomType()
    flat = ogr.GT_Flatten(gtype)
    print("input:  %s" % path)
    print("  geometry: %s (%d features)" % (ogr.GeometryTypeToName(gtype), lyr.GetFeatureCount()))
    if flat not in ACCEPTED:
        sys.exit("error: FARSITE's ignition reader handles Point, MultiPoint, "
                 "PolyLine and Polygon; %s is not one of them"
                 % ogr.GeometryTypeToName(gtype))
    if ogr.GT_HasZ(gtype) or ogr.GT_HasM(gtype):
        print("  note: Z/M values present; FARSITE reads X/Y only, dropping them")

    src_srs = lyr.GetSpatialRef()
    rings, npoints, nlines = [], 0, 0
    for feat in lyr:
        geom = feat.GetGeometryRef()
        if geom is None:
            continue
        for g in flatten_geom(geom):
            gt = ogr.GT_Flatten(g.GetGeometryType())
            if gt == ogr.wkbPolygon:
                outer = g.GetGeometryRef(0)
                pts = [(outer.GetX(i), outer.GetY(i)) for i in range(outer.GetPointCount())]
                if g.GetGeometryCount() > 1:
                    print("  note: feature has %d interior ring(s); FARSITE's reader "
                          "treats every ring as a separate outward-burning fire, so "
                          "holes are dropped" % (g.GetGeometryCount() - 1))
                rings.append(pts)
            elif gt == ogr.wkbLineString:
                pts = [(g.GetX(i), g.GetY(i)) for i in range(g.GetPointCount())]
                nlines += 1
                rings.append(pts)
            elif gt == ogr.wkbPoint:
                npoints += 1
                if args.point_to_circle:
                    rings.append(make_circle(g.GetX(), g.GetY(), args.radius,
                                             args.vertices, ccw=True))
                else:
                    rings.append([(g.GetX(), g.GetY())])
    ds = None

    if nlines:
        print("  note: %d polyline feature(s) -- FARSITE treats these as line\n"
              "        sources, spreading outward from the line rather than from an\n"
              "        area. This needs the fsxwignt.cpp:191 fix (`k` -> `j`); on an\n"
              "        unfixed build the simulation never terminates. See\n"
              "        tests/test_line_ignition.sh." % nlines)

    if npoints and not args.point_to_circle:
        print("  WARNING: %d point ignition(s). FARSITE expands each to a 10-vertex\n"
              "           circle of radius 1 map unit (startsize is hard-coded), which\n"
              "           is far below the %g m cell size. Pass --point-to-circle to\n"
              "           emit %g m polygon circles instead." % (npoints, lcp.resx, args.radius))
    if not rings:
        sys.exit("error: no usable geometry found in %s" % path)
    return rings, src_srs, flat


def flatten_geom(geom):
    """Yield single-part geometries from a possibly multi-part one."""
    gt = ogr.GT_Flatten(geom.GetGeometryType())
    if gt in (ogr.wkbMultiPoint, ogr.wkbMultiLineString, ogr.wkbMultiPolygon,
              ogr.wkbGeometryCollection):
        for i in range(geom.GetGeometryCount()):
            for g in flatten_geom(geom.GetGeometryRef(i)):
                yield g
    else:
        yield geom


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--lcp", required=True, help="landscape file to place/check against")
    p.add_argument("-o", "--output", required=True, help="output .shp path")
    p.add_argument("--mode", choices=["circle", "check"], default="circle")
    p.add_argument("--center", nargs=2, type=float, metavar=("X", "Y"),
                   help="circle centre in landscape coordinates (default: landscape centre)")
    p.add_argument("--radius", type=float, default=200.0,
                   help="circle radius in map units (default: 200)")
    p.add_argument("--vertices", type=int, default=32,
                   help="vertices in the circle (default: 32)")
    p.add_argument("--winding", choices=["ccw", "cw"], default="ccw",
                   help="ring direction. ccw matches what FARSITE's arp() wants "
                        "natively; cw matches the ESRI shapefile spec. FARSITE "
                        "normalises either way (default: ccw)")
    p.add_argument("--input", help="existing shapefile to validate (mode=check)")
    p.add_argument("--point-to-circle", action="store_true",
                   help="convert point ignitions to --radius polygon circles")
    p.add_argument("--min-edge-cells", type=float, default=2.0,
                   help="require this much clearance from the landscape edge, in "
                        "cells (default: 2)")
    p.add_argument("--crs", help="CRS for the output .prj (e.g. EPSG:5070). "
                                 "Default: taken from crop_info.json or a sibling .prj")
    p.add_argument("--force", action="store_true",
                   help="write the output even if validation found problems")
    args = p.parse_args()

    lcp = Lcp(args.lcp)
    lcp.describe()

    srs, srs_src = find_crs(args.lcp, args.crs)
    if srs is not None:
        print("  CRS for output: %s (from %s)"
              % (srs.GetName() or "unnamed", srs_src))
    else:
        print("  CRS for output: none found; writing without a .prj. The .lcp format "
              "stores no CRS, so pass --crs if you want one.")

    if args.mode == "circle":
        if args.radius <= 0:
            sys.exit("error: --radius must be positive")
        if args.vertices < 3:
            sys.exit("error: --vertices must be at least 3")
        cx, cy = args.center if args.center else lcp.center
        if not args.center:
            print("\nusing the landscape centre as the ignition centre")
        rings = [make_circle(cx, cy, args.radius, args.vertices,
                             ccw=(args.winding == "ccw"))]
        print("circle: centre (%.1f, %.1f), radius %g, %d vertices, %s"
              % (cx, cy, args.radius, args.vertices, args.winding))
        out_srs = srs
        out_type = ogr.wkbPolygon
    else:
        if not args.input:
            p.error("--mode check needs --input")
        rings, src_srs, out_type = read_and_clean(args.input, lcp, args)
        out_srs = src_srs or srs
        if src_srs is not None and srs is not None and not src_srs.IsSame(srs):
            print("  WARNING: the shapefile's CRS differs from the landscape's. "
                  "FARSITE does no reprojection at all -- coordinates are compared "
                  "raw -- so they must already match.")

    # Normalise: close rings, drop duplicate vertices, fix winding.
    cleaned = []
    for pts in rings:
        if is_closed(pts) or (args.mode == "circle"):
            pts = dedupe(close_ring(pts))
            if len(pts) < 4:
                print("  note: dropping a ring that collapsed to < 3 distinct vertices")
                continue
            want_ccw = args.winding == "ccw"
            if (signed_area(pts) > 0) != want_ccw:
                pts = pts[::-1]
            cleaned.append(pts)
        else:
            # Open line or point: order is meaningful, leave it alone.
            cleaned.append(dedupe(pts) if len(pts) > 1 else pts)
    rings = cleaned
    if not rings:
        sys.exit("error: nothing left after cleaning")

    for i, r in enumerate(rings):
        if is_closed(r):
            a = abs(signed_area(r))
            print("  ring %d: %d vertices, area %.0f m2 (%.2f ha), perimeter %.0f m, %s"
                  % (i, len(r) - 1, a, a / 10000.0, perimeter(r),
                     "ccw" if signed_area(r) > 0 else "cw"))
        elif len(r) > 1:
            print("  line %d: %d vertices, length %.0f m" % (i, len(r), perimeter(r)))
        else:
            print("  point %d: (%.1f, %.1f)" % (i, r[0][0], r[0][1]))

    problems, warnings = report_footprint(lcp, rings)

    margin_needed = args.min_edge_cells * lcp.resx
    xs = [pt[0] for r in rings for pt in r]
    ys = [pt[1] for r in rings for pt in r]
    corner_margin = min(lcp.edge_margin(x, y) for x, y in
                        [(min(xs), min(ys)), (min(xs), max(ys)),
                         (max(xs), min(ys)), (max(xs), max(ys))])
    if 0 <= corner_margin < margin_needed:
        warnings.append("only %.1f m (%.1f cells) of clearance to the landscape edge; "
                        "--min-edge-cells asks for %.1f. A fire that reaches the edge "
                        "simply stops, which truncates the run."
                        % (corner_margin, corner_margin / lcp.resx, args.min_edge_cells))

    for w in warnings:
        print("\n  WARNING: %s" % w)
    for e in problems:
        print("\n  PROBLEM: %s" % e)

    if problems and not args.force:
        sys.exit("\nrefusing to write %s. Fix the above, or pass --force."
                 % args.output)

    if all(is_closed(r) for r in rings):
        write_polygons(args.output, rings, out_srs)
    elif any(len(r) >= 2 for r in rings):
        # Open geometry: keep it a PolyLine so FARSITE treats it as a line source.
        write_lines(args.output, [r for r in rings if len(r) >= 2], out_srs)
        print("\nkept as PolyLine: FARSITE reads this as a line source, not an area.")
    else:
        sys.exit("error: only single points left after cleaning. FARSITE expands a "
                 "point to a 1 map-unit circle, which will not ignite a %g m cell; "
                 "re-run with --point-to-circle." % lcp.resx)

    print("\nwrote %s" % args.output)
    for ext in (".shx", ".dbf", ".prj"):
        side = os.path.splitext(args.output)[0] + ext
        if os.path.exists(side):
            print("      %s" % os.path.basename(side))
    print("\nUse it as the third token of the FARSITE command-file line:")
    print("  <lcp> <input> %s <barrier|0> <outbase> <type>" % os.path.abspath(args.output))
    return 0


if __name__ == "__main__":
    sys.exit(main())
