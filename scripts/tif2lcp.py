#!/usr/bin/env python3
"""Build a FARSITE landscape (.lcp) file from GeoTIFF inputs.

LANDFIRE stopped distributing .lcp in 2024 and now ships landscapes as
(multi-band) GeoTIFFs, but FARSITE only reads .lcp. This converts them.

Why not just `gdal_translate -of LCP`? GDAL's LCP writer produces correct pixel
data but a broken header statistics block: it writes each theme's category list
with a +32768 bias (it cannot re-read its own output) and folds the -9999 nodata
fill into the per-theme minima. That matters -- a wrong fuel-model category list
makes FARSITE's conditioning code build no fuel-moisture keys at all, so every
conditioned dead fuel moisture silently comes back 0.0 and the whole landscape
burns as if bone dry. This writes the header itself, field by field.

Usage
-----
  # LANDFIRE-style multi-band stack (bands in canonical LCP order)
  tif2lcp.py --stack landscape.tif -o out.lcp --latitude 45

  # or one GeoTIFF per theme
  tif2lcp.py --elevation e.tif --slope s.tif --aspect a.tif \
             --fuel f.tif --cover c.tif -o out.lcp --latitude 45

Band order is the LCP order: elevation, slope, aspect, fuel model, canopy cover,
then optionally canopy height, canopy base height, canopy bulk density, then
optionally duff and coarse woody debris. The five basic themes are required;
crown fuels (3) and ground fuels (2) are all-or-nothing groups.

Units default to what LANDFIRE ships. Override them if your rasters differ --
they are recorded in the header, not applied to the data, so a wrong unit code
silently misinterprets every value. See --help for the options and
`doc/lcp_format.md` for the semantics.
"""

import argparse
import os
import struct
import sys

import contextlib
import io

try:
    from osgeo import gdal
except ImportError:
    sys.exit("error: GDAL's Python bindings are required (pip install gdal, or conda install gdal)")

# UseExceptions() probes for the optional gdal_array module, which prints a bare
# "ModuleNotFoundError: No module named 'numpy'" to stderr when numpy is absent.
# This script deliberately avoids numpy (see read_theme), so swallow that probe
# rather than leaving what looks like a fatal error at the top of every run.
with contextlib.redirect_stderr(io.StringIO()):
    gdal.UseExceptions()

NODATA = -9999
HEADER_SIZE = 7316
NUM_CATS = 100

# theme key, CLI flag, human name, header unit field, default unit code
THEMES = [
    ("elevation", "elevation", "Elevation",          "EUnits",   0),
    ("slope",     "slope",     "Slope",              "SUnits",   0),
    ("aspect",    "aspect",    "Aspect",             "AUnits",   2),
    ("fuel",      "fuel",      "Fuel models",        "FOptions", 0),
    ("cover",     "cover",     "Canopy cover",       "CUnits",   1),
    ("height",    "height",    "Canopy height",      "HUnits",   3),
    ("base",      "base",      "Canopy base height", "BUnits",   3),
    ("density",   "density",   "Canopy bulk density","PUnits",   3),
    ("duff",      "duff",      "Duff",               "DUnits",   0),
    ("woody",     "woody",     "Coarse woody debris","WOptions", 0),
]
BASIC = ["elevation", "slope", "aspect", "fuel", "cover"]
CROWN = ["height", "base", "density"]
GROUND = ["duff", "woody"]

# Themes whose values are class codes: resampling must never interpolate them.
CATEGORICAL = {"fuel", "aspect", "cover", "woody"}

UNIT_CHOICES = {
    "EUnits":   {"meters": 0, "feet": 1},
    "SUnits":   {"degrees": 0, "percent": 1},
    "AUnits":   {"grass-categories": 0, "grass-degrees": 1, "azimuth-degrees": 2},
    "CUnits":   {"categories": 0, "percent": 1},
    "HUnits":   {"meters": 1, "feet": 2, "meters-x10": 3, "feet-x10": 4},
    "BUnits":   {"meters": 1, "feet": 2, "meters-x10": 3, "feet-x10": 4},
    "PUnits":   {"kg-per-m3": 1, "lb-per-ft3": 2, "kg-per-m3-x100": 3, "lb-per-ft3-x1000": 4},
    "DUnits":   {"class": 0, "tons-per-acre-x10": 1, "mg-per-ha-x10": 2},
    "FOptions": {"no-custom-no-file": 0, "custom-only": 1, "file-only": 2, "custom-and-file": 3},
    "WOptions": {"class": 0},
}


def theme_stats(values, count):
    """LCP per-theme statistics: (lo, hi, num, cats[100]).

    The category list is 1-based: slot 0 is unused, the distinct values occupy
    slots 1..num sorted ascending, and num is -1 when there are more than 99
    distinct values. lo/hi ignore nodata. This mirrors what FARSITE's own
    LandscapeTheme::FillCats/SortCats produce, which is the de facto spec.
    """
    valid = [v for v in values if v != NODATA and v >= 0]
    lo = min(valid) if valid else 0
    hi = max(values) if values else 0

    distinct = sorted(set(values))
    cats = [0] * NUM_CATS
    if len(distinct) > NUM_CATS - 2:
        num = -1
    else:
        num = len(distinct)
        for i, v in enumerate(distinct):
            cats[i + 1] = v
    return lo, hi, num, cats


def open_theme_bands(args):
    """Return {theme: (gdal.Band, dataset)} for every theme the user supplied."""
    bands = {}
    if args.stack:
        ds = gdal.Open(args.stack)
        n = ds.RasterCount
        order = BASIC + (CROWN if n >= 8 else []) + (GROUND if n in (7, 10) else [])
        if n == 7:
            order = BASIC + GROUND
        if n not in (5, 7, 8, 10):
            sys.exit("error: --stack has %d bands; expected 5, 7, 8 or 10 "
                     "(see band order in --help)" % n)
        for i, key in enumerate(order):
            bands[key] = (ds.GetRasterBand(i + 1), ds)
    for key, flag, _, _, _ in THEMES:
        path = getattr(args, flag)
        if path:
            if key in bands:
                sys.exit("error: theme '%s' given both in --stack and as --%s" % (key, flag))
            ds = gdal.Open(path)
            if ds.RasterCount != 1:
                sys.exit("error: --%s (%s) has %d bands; expected 1"
                         % (flag, path, ds.RasterCount))
            bands[key] = (ds.GetRasterBand(1), ds)
    return bands


def check_grids(bands):
    """LCP stores one grid for all themes, so every input must already agree."""
    ref_key = BASIC[0]
    rb, rds = bands[ref_key]
    ref = (rds.RasterXSize, rds.RasterYSize, rds.GetGeoTransform())
    for key, (b, ds) in bands.items():
        got = (ds.RasterXSize, ds.RasterYSize, ds.GetGeoTransform())
        if got[:2] != ref[:2] or any(abs(a - c) > 1e-6 for a, c in zip(got[2], ref[2])):
            sys.exit("error: theme '%s' is not on the same grid as '%s'.\n"
                     "  %s: %dx%d gt=%s\n  %s: %dx%d gt=%s\n"
                     "LCP requires one shared extent/resolution/alignment. Align them first, e.g.\n"
                     "  gdalwarp -te <xmin ymin xmax ymax> -tr <xres> <yres> -r near in.tif out.tif"
                     % (key, ref_key, key, got[0], got[1], got[2],
                        ref_key, ref[0], ref[1], ref[2]))
    gt = ref[2]
    if abs(gt[2]) > 1e-9 or abs(gt[4]) > 1e-9:
        sys.exit("error: rotated/skewed geotransform; LCP requires a north-up grid")
    if gt[5] > 0:
        sys.exit("error: south-up raster; LCP requires north-up (negative y pixel size)")
    return rds.RasterXSize, rds.RasterYSize, gt


def read_theme(band, key, ncols, nrows):
    """Read a band as a flat list of ints, normalising nodata to -9999.

    Uses ReadRaster rather than ReadAsArray so the script needs only GDAL, not
    numpy, and reads a row at a time to keep memory bounded on large
    landscapes. Values come back as float64 so that scaled or floating-point
    source rasters round correctly.
    """
    nd = band.GetNoDataValue()
    unpack = struct.Struct("<%dd" % ncols).unpack
    vals = []
    for row in range(nrows):
        raw = band.ReadRaster(0, row, ncols, 1, ncols, 1, gdal.GDT_Float64)
        for v in unpack(raw):
            if nd is not None and v == nd:
                iv = NODATA
            else:
                iv = int(round(v))
            if iv < -32768 or iv > 32767:
                sys.exit("error: theme '%s' has value %d at row %d, outside the "
                         "int16 range LCP stores; check units/scaling"
                         % (key, iv, row))
            vals.append(iv)
    return vals


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-o", "--output", required=True, help="output .lcp path")
    p.add_argument("--stack", help="multi-band GeoTIFF with bands in LCP order")
    for key, flag, name, _, _ in THEMES:
        p.add_argument("--%s" % flag, help="single-band GeoTIFF for %s" % name.lower())
    p.add_argument("--latitude", type=int,
                   help="latitude in whole degrees; drives solar radiation. "
                        "Taken from the raster CRS when omitted and derivable.")
    p.add_argument("--description", default="", help="free-text header description (max 511 chars)")
    p.add_argument("--linear-unit", choices=["meters", "feet"], default="meters",
                   help="horizontal unit of the grid (default: meters)")
    for field, choices in UNIT_CHOICES.items():
        dflt = next(d for _, _, _, f, d in THEMES if f == field)
        dname = next(k for k, v in choices.items() if v == dflt)
        p.add_argument("--%s" % field.lower(), dest=field, choices=sorted(choices),
                       default=dname, help="units for %s (default: %s)" % (field, dname))
    args = p.parse_args()

    if not args.stack and not any(getattr(args, f) for _, f, _, _, _ in THEMES):
        p.error("give either --stack or per-theme GeoTIFFs")

    bands = open_theme_bands(args)

    missing = [k for k in BASIC if k not in bands]
    if missing:
        sys.exit("error: missing required theme(s): %s" % ", ".join(missing))
    have_crown = [k for k in CROWN if k in bands]
    if have_crown and len(have_crown) != 3:
        sys.exit("error: crown fuels are all-or-nothing; got %s, need all of %s"
                 % (", ".join(have_crown), ", ".join(CROWN)))
    have_ground = [k for k in GROUND if k in bands]
    if have_ground and len(have_ground) != 2:
        sys.exit("error: ground fuels are all-or-nothing; got %s, need both of %s"
                 % (", ".join(have_ground), ", ".join(GROUND)))

    ncols, nrows, gt = check_grids(bands)
    crown_fuels = 21 if have_crown else 20
    ground_fuels = 21 if have_ground else 20
    order = BASIC + (CROWN if have_crown else []) + (GROUND if have_ground else [])
    numvals = len(order)

    latitude = args.latitude
    if latitude is None:
        rds = bands[BASIC[0]][1]
        try:
            from osgeo import osr
            srs = osr.SpatialReference(wkt=rds.GetProjection())
            if srs.IsProjected() or srs.IsGeographic():
                tgt = osr.SpatialReference(); tgt.ImportFromEPSG(4326)
                tgt.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
                tr = osr.CoordinateTransformation(srs, tgt)
                cx = gt[0] + gt[1] * ncols / 2.0
                cy = gt[3] + gt[5] * nrows / 2.0
                latitude = int(round(tr.TransformPoint(cx, cy)[1]))
        except Exception:
            latitude = None
    if latitude is None:
        sys.exit("error: --latitude is required (the raster has no usable CRS to derive it from).\n"
                 "FARSITE uses it for solar radiation, so it must be right.")
    if not -90 <= latitude <= 90:
        sys.exit("error: --latitude %d out of range" % latitude)

    print("grid:   %d cols x %d rows, cell %.6g x %.6g" % (ncols, nrows, gt[1], abs(gt[5])))
    print("themes: %s  (NumVals=%d)" % (", ".join(order), numvals))
    print("lat:    %d" % latitude)

    # Read every theme and compute its statistics.
    data, stats = {}, {}
    for key in order:
        band, _ = bands[key]
        vals = read_theme(band, key, ncols, nrows)
        if len(vals) != ncols * nrows:
            sys.exit("error: theme '%s' returned %d cells, expected %d"
                     % (key, len(vals), ncols * nrows))
        data[key] = vals
        stats[key] = theme_stats(vals, ncols * nrows)
        lo, hi, num, _ = stats[key]
        print("  %-9s lo=%-8d hi=%-8d classes=%s" % (key, lo, hi,
              num if num >= 0 else ">99 (unclassified)"))

    west, north = gt[0], gt[3]
    east = west + gt[1] * ncols
    south = north + gt[5] * nrows

    # ---- header, written field by field in on-disk order ----
    # Note sizeof(headdata) in memory is 7328 because of double alignment
    # padding; the on-disk layout has none, so it cannot be packed as one struct.
    h = bytearray()
    h += struct.pack("<3i", crown_fuels, ground_fuels, latitude)
    h += struct.pack("<4d", west, east, south, north)   # loeast hieast lonorth hinorth
    for key, _, _, _, _ in THEMES:
        if key in stats:
            lo, hi, num, cats = stats[key]
        else:
            lo, hi, num, cats = 0, 0, 0, [0] * NUM_CATS
        h += struct.pack("<3i", lo, hi, num)
        h += struct.pack("<%di" % NUM_CATS, *cats)
    h += struct.pack("<2i", ncols, nrows)
    h += struct.pack("<4d", east, west, north, south)
    h += struct.pack("<i", 0 if args.linear_unit == "meters" else 1)
    h += struct.pack("<2d", abs(gt[1]), abs(gt[5]))
    for key, _, _, field, _ in THEMES:
        h += struct.pack("<h", UNIT_CHOICES[field][getattr(args, field)])
    for key, flag, _, _, _ in THEMES:
        src = getattr(args, flag) or (args.stack if key in bands else "")
        h += os.path.basename(src or "").encode()[:255].ljust(256, b"\0")
    h += args.description.encode()[:511].ljust(512, b"\0")
    assert len(h) == HEADER_SIZE, "header is %d bytes, expected %d" % (len(h), HEADER_SIZE)

    # ---- interleaved cell records, row-major from the north edge down ----
    planes = [data[k] for k in order]
    body = bytearray()
    pack = struct.Struct("<%dh" % numvals).pack
    for i in range(ncols * nrows):
        body += pack(*[pl[i] for pl in planes])

    with open(args.output, "wb") as f:
        f.write(h)
        f.write(body)

    size = HEADER_SIZE + ncols * nrows * numvals * 2
    print("wrote %s (%d bytes)" % (args.output, size))
    print("\nVerify with:  gdalinfo %s" % args.output)
    print("If FARSITE rejects it or results look wrong, re-run FARSITE with")
    print("--recompute-lcp-stats to rebuild the statistics block from the raster.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
