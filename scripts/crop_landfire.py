#!/usr/bin/env python3
"""Crop LANDFIRE CONUS GeoTIFFs to a square window around a fire centroid.

Produces one single-band GeoTIFF per LCP theme, all on an identical grid, ready
to hand to scripts/tif2lcp.py. Themes that were never downloaded are synthesised
as constant rasters on that same grid so the set is always complete.

Inputs are two CSVs.

Availability manifest -- which CONUS mosaics are on disk:

    survey_year,band,gcs_path
    2025,fuel,/path/LF2025_FBFM40_CONUS.tif
    2020,elevation,/path/LF2020_Elev_CONUS.tif

Training manifest -- the fires, with centroids in EPSG:4326:

    id,fire_name,date,data_type,status,error_msg,centroid_x,centroid_y
    0,Dish,2020-08-20,polygon,ready,,-122.17931835079123,37.408402398277026

Example
-------
    # 2 km box (67 * 30 m) at native resolution, one fire
    python3 scripts/crop_landfire.py \
        --available-manifest available.csv \
        --training-manifest training.csv \
        --fire-id 0 --size-cells 67 -o cases/

    # same footprint resampled to 90 m
    python3 scripts/crop_landfire.py ... --size-cells 67 --resolution 90

The window side is given in multiples of the native LANDFIRE cell (30 m), so
--size-cells 67 is 2010 m on a side. The window is snapped to the source raster
grid, so at native resolution the crop is a pure subset with no resampling.
"""

import argparse
import csv
import json
import os
import re
import shlex
import struct
import sys

try:
    from osgeo import gdal, osr
except ImportError:
    sys.exit("error: GDAL's Python bindings are required (pip install gdal, or conda install gdal)")

# UseExceptions() probes the optional gdal_array module, which prints a bare
# "ModuleNotFoundError: No module named 'numpy'" when numpy is absent. This
# script deliberately avoids numpy, so swallow that probe.
import contextlib
import io
with contextlib.redirect_stderr(io.StringIO()):
    gdal.UseExceptions()
    osr.UseExceptions()

class WindowError(Exception):
    """A fire cannot be cropped; skip it rather than aborting the batch."""


NODATA = -9999
NATIVE_RES = 30.0          # LANDFIRE CONUS cell size, metres

# LCP theme order. Each entry carries the default value used when the theme was
# not downloaded, the resampling algorithm, and the tif2lcp unit flag implied by
# the LANDFIRE product that normally supplies it.
#
# 'near' is mandatory for fuel (class codes: interpolating invents models that
# do not exist) and for aspect (azimuth degrees wrap at 360, so averaging 359
# and 1 gives 180 -- due south instead of due north).
THEMES = {
    "elevation": dict(default=0,   alg="bilinear", unit=("eunits", "meters")),
    "slope":     dict(default=0,   alg="bilinear", unit=("sunits", "degrees")),
    "aspect":    dict(default=0,   alg="near",     unit=("aunits", "azimuth-degrees")),
    "fuel":      dict(default=101, alg="near",     unit=("foptions", "no-custom-no-file")),
    "cover":     dict(default=0,   alg="bilinear", unit=("cunits", "percent")),
    "height":    dict(default=0,   alg="bilinear", unit=("hunits", "meters-x10")),
    "base":      dict(default=0,   alg="bilinear", unit=("bunits", "meters-x10")),
    "density":   dict(default=0,   alg="bilinear", unit=("punits", "kg-per-m3-x100")),
    "duff":      dict(default=0,   alg="near",     unit=("dunits", "class")),
    "woody":     dict(default=0,   alg="near",     unit=("woptions", "class")),
}
BASIC = ["elevation", "slope", "aspect", "fuel", "cover"]
CROWN = ["height", "base", "density"]
GROUND = ["duff", "woody"]

CREATE_OPTS = ["COMPRESS=DEFLATE", "TILED=YES"]


def read_csv(path):
    """Read a CSV, stripping whitespace from headers and every value.

    Both manifests in practice have trailing spaces on the header line and
    padding after the paths, which would otherwise become part of a filename.
    """
    with open(path, newline="") as fh:
        rows = []
        for row in csv.DictReader(fh):
            rows.append({(k or "").strip(): (v or "").strip() for k, v in row.items()})
        return rows


def load_available(path, prefer_year=None):
    """{band: {path, year, all_years}} picking one raster per band."""
    by_band = {}
    for r in read_csv(path):
        band, src = r.get("band", ""), r.get("gcs_path", "")
        if not band or not src:
            continue
        if band not in THEMES:
            print("  warning: ignoring unknown band '%s' in availability manifest" % band)
            continue
        try:
            year = int(r.get("survey_year") or 0)
        except ValueError:
            year = 0
        by_band.setdefault(band, []).append((year, src))

    chosen = {}
    for band, entries in by_band.items():
        entries.sort(key=lambda e: e[0], reverse=True)
        pick = None
        if prefer_year is not None:
            pick = next((e for e in entries if e[0] == prefer_year), None)
            if pick is None:
                print("  warning: %s has no %d survey; using %d"
                      % (band, prefer_year, entries[0][0]))
        pick = pick or entries[0]
        missing = [p for _, p in entries if not os.path.exists(p)]
        if not os.path.exists(pick[1]):
            print("  warning: %s listed at %s but the file is not readable; "
                  "treating as missing" % (band, pick[1]))
            continue
        chosen[band] = dict(path=pick[1], year=pick[0],
                            all_years=sorted({y for y, _ in entries}, reverse=True))
        if missing:
            print("  note: %d other %s entr%s unreadable"
                  % (len(missing), band, "y" if len(missing) == 1 else "ies"))
    return chosen


def detect_units(band, path):
    """Override the assumed unit flag when the filename says otherwise.

    LANDFIRE ships slope as both SlpD (degrees) and SlpP (percent) under
    otherwise identical naming, and the unit code is recorded in the LCP header
    rather than applied to the data -- so getting it wrong silently
    misinterprets every cell.
    """
    name = os.path.basename(path)
    if band == "slope":
        if re.search(r"SlpP", name, re.I):
            return ("sunits", "percent")
        if re.search(r"SlpD", name, re.I):
            return ("sunits", "degrees")
        print("  warning: cannot tell slope units from '%s'; assuming degrees. "
              "Pass --sunits to tif2lcp if that is wrong." % name)
    return THEMES[band]["unit"]


def slug(text):
    s = re.sub(r"[^A-Za-z0-9]+", "_", text or "").strip("_")
    return s or "unnamed"


def reference_grid(available):
    """Pick a raster to define the target CRS and grid origin."""
    for band in BASIC + CROWN + GROUND:
        if band in available:
            ds = gdal.Open(available[band]["path"])
            return band, ds
    return None, None


def compute_window(ds, lon, lat, size_cells, out_res, native_res):
    """Square window around (lon, lat), snapped to the source raster grid.

    Returns (xmin, ymin, xmax, ymax, ncols, nrows, wkt).
    """
    wkt = ds.GetProjection()
    if not wkt:
        raise WindowError("reference raster has no CRS; cannot place the centroid")
    dst = osr.SpatialReference(wkt=wkt)
    src = osr.SpatialReference()
    src.ImportFromEPSG(4326)
    src.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    dst.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    cx, cy, _ = osr.CoordinateTransformation(src, dst).TransformPoint(lon, lat)

    side = size_cells * native_res
    ncols = nrows = int(round(side / out_res))
    if ncols < 1:
        raise WindowError("window is %.1f m at %.1f m resolution, which is less "
                          "than one cell" % (side, out_res))

    # Snap the top-left corner to the source grid so that, when out_res equals
    # the source resolution, the crop is an exact pixel subset (no resampling).
    gt = ds.GetGeoTransform()
    ox, oy = gt[0], gt[3]
    xmin = ox + round((cx - side / 2.0 - ox) / out_res) * out_res
    ymax = oy - round((oy - (cy + side / 2.0)) / out_res) * out_res
    xmax = xmin + ncols * out_res
    ymin = ymax - nrows * out_res

    # The centroid must actually land inside the mosaic.
    rx0, ry1 = gt[0], gt[3]
    rx1 = rx0 + gt[1] * ds.RasterXSize
    ry0 = ry1 + gt[5] * ds.RasterYSize
    if not (rx0 <= cx <= rx1 and ry0 <= cy <= ry1):
        raise WindowError("centroid (%.6f, %.6f) projects to (%.1f, %.1f), outside "
                          "the reference raster extent x[%.1f, %.1f] y[%.1f, %.1f]"
                          % (lon, lat, cx, cy, rx0, rx1, ry0, ry1))
    if not (rx0 <= xmin and xmax <= rx1 and ry0 <= ymin and ymax <= ry1):
        print("  warning: the window extends past the edge of the CONUS mosaic; "
              "outside cells will be nodata (-9999)")
    return xmin, ymin, xmax, ymax, ncols, nrows, wkt, (cx, cy)


def warp_theme(src_path, out_path, bounds, out_res, wkt, alg, want_cols, want_rows):
    xmin, ymin, xmax, ymax = bounds
    gdal.Warp(out_path, src_path, format="GTiff",
              outputBounds=(xmin, ymin, xmax, ymax), outputBoundsSRS=wkt,
              dstSRS=wkt, xRes=out_res, yRes=out_res,
              resampleAlg=alg, dstNodata=NODATA,
              outputType=gdal.GDT_Int16, creationOptions=CREATE_OPTS)
    ds = gdal.Open(out_path)
    if (ds.RasterXSize, ds.RasterYSize) != (want_cols, want_rows):
        sys.exit("error: %s came out %dx%d, expected %dx%d"
                 % (out_path, ds.RasterXSize, ds.RasterYSize, want_cols, want_rows))


def write_constant(out_path, value, bounds, out_res, wkt, ncols, nrows):
    """Synthesise a constant Int16 raster on the target grid (no numpy)."""
    xmin, ymax = bounds[0], bounds[3]
    drv = gdal.GetDriverByName("GTiff")
    ds = drv.Create(out_path, ncols, nrows, 1, gdal.GDT_Int16, options=CREATE_OPTS)
    ds.SetGeoTransform((xmin, out_res, 0.0, ymax, 0.0, -out_res))
    ds.SetProjection(wkt)
    band = ds.GetRasterBand(1)
    band.SetNoDataValue(NODATA)
    row = struct.pack("=%dh" % ncols, *([value] * ncols))
    for r in range(nrows):
        band.WriteRaster(0, r, ncols, 1, row, buf_type=gdal.GDT_Int16)
    band.FlushCache()
    ds = None


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--available-manifest", required=True,
                   help="CSV of downloaded CONUS mosaics (survey_year,band,gcs_path)")
    p.add_argument("--training-manifest", required=True,
                   help="CSV of fires (id,fire_name,...,centroid_x,centroid_y in EPSG:4326)")
    p.add_argument("-o", "--outdir", required=True, help="output directory")
    p.add_argument("--size-cells", type=int, required=True,
                   help="window side length in multiples of the native cell "
                        "(%g m), e.g. 67 -> %g m" % (NATIVE_RES, 67 * NATIVE_RES))
    p.add_argument("--resolution", type=float,
                   help="output cell size in metres (default: native, %g)" % NATIVE_RES)
    p.add_argument("--native-res", type=float, default=NATIVE_RES,
                   help="native LANDFIRE cell size (default: %g)" % NATIVE_RES)
    p.add_argument("--fire-id", action="append", default=[],
                   help="fire id to process; repeatable. Default: all matching --status")
    p.add_argument("--fire-name", action="append", default=[],
                   help="fire name to process; repeatable")
    p.add_argument("--status", default="ready",
                   help="only process rows with this status ('any' to disable; default: ready)")
    p.add_argument("--prefer-year", type=int,
                   help="prefer this LANDFIRE survey year when several are available "
                        "(default: newest per band)")
    p.add_argument("--include-ground-fuels", action="store_true",
                   help="also emit duff and coarse woody debris. LANDFIRE does not "
                        "ship these, so they would be synthetic; omitting them makes "
                        "the LCP declare GroundFuels=absent, which is more honest.")
    p.add_argument("--resample-alg",
                   help="override the resampling algorithm for every theme. Per-theme "
                        "defaults are nearest for fuel and aspect, bilinear otherwise.")
    for band, spec in THEMES.items():
        p.add_argument("--default-%s" % band, type=int, default=spec["default"],
                       help="value for a synthesised %s raster (default: %d)"
                            % (band, spec["default"]))
    p.add_argument("--dry-run", action="store_true", help="report the plan and stop")
    args = p.parse_args()

    out_res = args.resolution or args.native_res
    wanted = BASIC + CROWN + (GROUND if args.include_ground_fuels else [])

    print("== availability ==")
    available = load_available(args.available_manifest, args.prefer_year)
    for band in wanted:
        if band in available:
            a = available[band]
            print("  %-9s %d  %s" % (band, a["year"], os.path.basename(a["path"])))
        else:
            print("  %-9s --    MISSING -> constant %d"
                  % (band, getattr(args, "default_%s" % band)))
    if not available:
        sys.exit("error: no usable rasters in the availability manifest")

    missing_basic = [b for b in BASIC if b not in available]
    if missing_basic:
        print("\n  WARNING: required theme(s) %s will be synthetic constants."
              % ", ".join(missing_basic))
        if "fuel" in missing_basic:
            print("  A uniform fuel model means the simulation tests plumbing, "
                  "not fire behaviour.")

    ref_band, ref_ds = reference_grid(available)
    print("\n  grid reference: %s" % ref_band)
    ref_gt = ref_ds.GetGeoTransform()
    print("  source cell: %g x %g" % (ref_gt[1], abs(ref_gt[5])))
    if abs(abs(ref_gt[1]) - args.native_res) > 1e-6:
        print("  warning: reference raster cell is %g m but --native-res is %g m"
              % (abs(ref_gt[1]), args.native_res))

    # Every source must share the reference CRS, or warp would silently
    # reproject and the "no resampling at native res" guarantee would be lost.
    ref_wkt = ref_ds.GetProjection()
    ref_srs = osr.SpatialReference(wkt=ref_wkt)
    for band, a in available.items():
        srs = osr.SpatialReference(wkt=gdal.Open(a["path"]).GetProjection())
        if not srs.IsSame(ref_srs):
            print("  note: %s is in a different CRS; it will be reprojected" % band)

    fires = read_csv(args.training_manifest)
    if args.status.lower() != "any":
        fires = [f for f in fires if f.get("status", "").lower() == args.status.lower()]
    if args.fire_id:
        want = {s.strip() for s in args.fire_id}
        fires = [f for f in fires if f.get("id", "") in want]
    if args.fire_name:
        want = {s.strip().lower() for s in args.fire_name}
        fires = [f for f in fires if f.get("fire_name", "").lower() in want]
    if not fires:
        sys.exit("error: no fires selected (check --fire-id/--fire-name/--status)")

    print("\n== %d fire(s), %d x %d cells at %g m (%g m across) ==" % (
        len(fires), int(round(args.size_cells * args.native_res / out_res)),
        int(round(args.size_cells * args.native_res / out_res)), out_res,
        args.size_cells * args.native_res))

    written, skipped = [], []
    for f in fires:
        fid, name = f.get("id", "?"), f.get("fire_name", "")
        try:
            lon, lat = float(f["centroid_x"]), float(f["centroid_y"])
        except (KeyError, ValueError):
            print("  SKIP id=%s %s: unreadable centroid" % (fid, name))
            continue
        if not (-180 <= lon <= 180 and -90 <= lat <= 90):
            print("  SKIP id=%s %s: centroid (%s, %s) is not valid EPSG:4326"
                  % (fid, name, lon, lat))
            continue

        case = "%s_%s" % (fid, slug(name))
        dest = os.path.join(args.outdir, case)
        print("\n-- %s  (%.5f, %.5f)" % (case, lon, lat))

        try:
            xmin, ymin, xmax, ymax, ncols, nrows, wkt, (cx, cy) = compute_window(
                ref_ds, lon, lat, args.size_cells, out_res, args.native_res)
        except WindowError as e:
            print("   SKIP: %s" % e)
            skipped.append((case, str(e)))
            continue
        bounds = (xmin, ymin, xmax, ymax)
        print("   centre in raster CRS: (%.1f, %.1f)" % (cx, cy))
        print("   window: x[%.1f, %.1f] y[%.1f, %.1f]  %dx%d" % (xmin, xmax, ymin, ymax, ncols, nrows))
        if args.dry_run:
            continue

        os.makedirs(dest, exist_ok=True)
        units, provenance = {}, {}
        for band in wanted:
            out_path = os.path.join(dest, "%s.tif" % band)
            if band in available:
                alg = args.resample_alg or THEMES[band]["alg"]
                if abs(out_res - abs(ref_gt[1])) < 1e-6:
                    alg = "near"   # 1:1 subset; nothing to interpolate
                warp_theme(available[band]["path"], out_path, bounds, out_res,
                           wkt, alg, ncols, nrows)
                key, val = detect_units(band, available[band]["path"])
                provenance[band] = dict(source=available[band]["path"],
                                        year=available[band]["year"], resample=alg)
                print("   %-9s cropped (%s)" % (band, alg))
            else:
                val0 = getattr(args, "default_%s" % band)
                write_constant(out_path, val0, bounds, out_res, wkt, ncols, nrows)
                key, val = THEMES[band]["unit"]
                provenance[band] = dict(source=None, constant=val0)
                print("   %-9s synthesised = %d" % (band, val0))
            units[key] = val

        # tif2lcp needs a latitude in whole degrees; the centroid is exact.
        lat_deg = int(round(lat))
        # Every argument is shell-quoted: LANDFIRE paths routinely contain
        # spaces (e.g. /Volumes/Extreme SSD/...) and the description does too.
        pairs = [("--%s" % band, os.path.join(dest, "%s.tif" % band)) for band in wanted]
        pairs.append(("-o", os.path.join(dest, "%s.lcp" % case)))
        pairs.append(("--latitude", str(lat_deg)))
        pairs += [("--%s" % k, v) for k, v in sorted(units.items())]
        pairs.append(("--description", "LANDFIRE crop %s" % case))

        lines = ["python3 " + shlex.quote(os.path.join("scripts", "tif2lcp.py"))]
        lines += ["  %s %s" % (flag, shlex.quote(val)) for flag, val in pairs]

        script_path = os.path.join(dest, "make_lcp.sh")
        with open(script_path, "w") as fh:
            fh.write("#!/bin/sh\n"
                     "# Generated by crop_landfire.py. Run from the repo root.\n"
                     "set -e\n")
            fh.write(" \\\n".join(lines) + "\n")
        os.chmod(script_path, 0o755)

        with open(os.path.join(dest, "crop_info.json"), "w") as fh:
            json.dump(dict(id=fid, fire_name=name, date=f.get("date"),
                           centroid_4326=[lon, lat], centre_projected=[cx, cy],
                           bounds=dict(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax),
                           ncols=ncols, nrows=nrows, resolution=out_res,
                           crs_wkt=wkt, themes=provenance,
                           tif2lcp_units=units, latitude=lat_deg), fh, indent=2)

        written.append(dest)
        print("   -> %s" % dest)
        print("      build the .lcp with: %s/make_lcp.sh" % dest)

    print()
    if written:
        print("wrote %d case(s) under %s" % (len(written), args.outdir))
    if skipped:
        print("skipped %d fire(s):" % len(skipped))
        for case, why in skipped:
            print("  %s: %s" % (case, why))
    if not written and not args.dry_run:
        print("no cases written")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
