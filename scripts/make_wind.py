#!/usr/bin/env python3
"""Generate gridded wind fields (.atm + ASCII grids) for FARSITE.

A RAWS or .wnd stream gives FARSITE one wind speed and direction for the whole
landscape at each time. To vary wind in space you need the gridded-wind path:
an .atm file listing, per time step, a wind-speed grid and a wind-direction
grid. This writes those.

Modes
-----
  constant  one uniform speed/direction over the whole grid (default)
  noisy     the same field with Gaussian noise on speed and direction,
            optionally spatially correlated
  hrrr      ingest HRRR forecast winds -- not implemented yet, see load_hrrr()

Examples
--------
  # 5 m/s from the west (blowing west-to-east), 6 hours, matching an LCP
  python3 scripts/make_wind.py --lcp examples/Dish/0_Dish/0_Dish.lcp \
      -o examples/Dish/0_Dish/wind --start "04 26 1200"

  # take the time window straight from the FARSITE input file
  python3 scripts/make_wind.py --lcp case.lcp -o wind --input-file case.input

  # 15-minute steps with noise
  python3 scripts/make_wind.py --lcp case.lcp -o wind --start "04 26 1200" \
      --interval-minutes 15 --mode noisy --noise-speed 1.5 --noise-direction 20

Units
-----
FARSITE's gridded-wind reader is picky, and the .atm units keyword decides how
it reads the speed grids (`CAtmWindGrid::Load`):

  METRIC   speed is km/h at 10 m. FARSITE divides by 1.15 to get a 20 ft wind,
           then converts to mph.
  ENGLISH  speed is mph at 20 ft, used as-is. This is the default when the .atm
           has no units keyword.

`--speed` is given in `--speed-units` (m/s by default, the usual way a 10 m wind
is quoted) and converted to whichever the output demands. METRIC output is the
default because it lets FARSITE apply the 10 m -> 20 ft reduction rather than
having this script guess at it.

Direction
---------
`--direction` is the compass azimuth the wind blows **from**, which is both the
meteorological convention and what FARSITE expects: 270 means a westerly, i.e.
air moving toward the east. That is the default. Use `--toward` if you would
rather name the direction the wind is heading.
"""

import argparse
import math
import os
import random
import re
import struct
import sys

NODATA = -9999.0

# LCP header field offsets (see doc/lcp_format.md). Read directly so this
# script needs no GDAL for the common case of matching an existing .lcp.
LCP_HEADER_SIZE = 7316
LCP_OFF_NUMEAST = 4164
LCP_OFF_EASTUTM = 4172

MPS_TO_KPH = 3.6
MPS_TO_MPH = 2.2369362920544
KPH_TO_MPS = 1.0 / MPS_TO_KPH
MPH_TO_MPS = 1.0 / MPS_TO_MPH

# FARSITE's own 10 m -> 20 ft reduction, applied to METRIC grids in
# CAtmWindGrid::Load. Reproduced here only to report the effective wind.
TEN_M_TO_TWENTY_FT = 1.0 / 1.15


class Grid(object):
    """Target grid for the wind rasters."""

    def __init__(self, ncols, nrows, xll, yll, cellsize):
        self.ncols, self.nrows = ncols, nrows
        self.xll, self.yll, self.cellsize = xll, yll, cellsize

    @property
    def xur(self):
        return self.xll + self.ncols * self.cellsize

    @property
    def yur(self):
        return self.yll + self.nrows * self.cellsize

    def header(self):
        return ("ncols\t%d\n"
                "nrows\t%d\n"
                "xllcorner\t%.6f\n"
                "yllcorner\t%.6f\n"
                "cellsize\t%.6f\n"
                "NODATA_value\t%.6f\n"
                % (self.ncols, self.nrows, self.xll, self.yll,
                   self.cellsize, NODATA))

    def __str__(self):
        return ("%d x %d cells at %g, x[%.1f, %.1f] y[%.1f, %.1f]"
                % (self.ncols, self.nrows, self.cellsize,
                   self.xll, self.xur, self.yll, self.yur))


def read_lcp_extent(path):
    """Extent and cell count of a .lcp, straight from its header (no GDAL)."""
    with open(path, "rb") as fh:
        head = fh.read(LCP_HEADER_SIZE)
    if len(head) < LCP_HEADER_SIZE:
        sys.exit("error: %s is too short to be an .lcp" % path)
    numeast, numnorth = struct.unpack_from("<2i", head, LCP_OFF_NUMEAST)
    east, west, north, south = struct.unpack_from("<4d", head, LCP_OFF_EASTUTM)
    if numeast <= 0 or numnorth <= 0:
        sys.exit("error: %s reports a %dx%d grid; not a readable .lcp"
                 % (path, numeast, numnorth))
    # FARSITE derives resolution from the corners, not from XResol/YResol.
    resx = (east - west) / numeast
    resy = (north - south) / numnorth
    return numeast, numnorth, west, south, east, north, resx, resy


def read_raster_extent(path):
    """Extent of any GDAL-readable raster, used for --like."""
    try:
        import contextlib
        import io
        from osgeo import gdal
        with contextlib.redirect_stderr(io.StringIO()):
            gdal.UseExceptions()
    except ImportError:
        sys.exit("error: --like needs GDAL's Python bindings; use --lcp instead "
                 "to read an .lcp header directly")
    ds = gdal.Open(path)
    gt = ds.GetGeoTransform()
    if abs(gt[2]) > 1e-9 or abs(gt[4]) > 1e-9:
        sys.exit("error: %s has a rotated geotransform" % path)
    west, north = gt[0], gt[3]
    east = west + gt[1] * ds.RasterXSize
    south = north + gt[5] * ds.RasterYSize
    return (ds.RasterXSize, ds.RasterYSize, west, south, east, north,
            abs(gt[1]), abs(gt[5]))


def build_grid(args):
    if args.lcp:
        ncols, nrows, west, south, east, north, resx, resy = read_lcp_extent(args.lcp)
        src = args.lcp
    else:
        ncols, nrows, west, south, east, north, resx, resy = read_raster_extent(args.like)
        src = args.like
    print("reference: %s" % src)
    print("  %d x %d cells, res %g x %g, x[%.1f, %.1f] y[%.1f, %.1f]"
          % (ncols, nrows, resx, resy, west, east, south, north))
    if abs(resx - resy) > 1e-6:
        print("  warning: non-square reference cells (%g x %g)" % (resx, resy))

    res = args.resolution or resx
    # FARSITE rejects a wind grid that does not cover the analysis area
    # (CWindGrids::CheckCoverage). Pad outward so edge rounding cannot bite.
    pad = args.pad_cells * res
    gxll = west - pad
    gyll = south - pad
    gncols = int(math.ceil((east + pad - gxll) / res))
    gnrows = int(math.ceil((north + pad - gyll) / res))
    grid = Grid(gncols, gnrows, gxll, gyll, res)
    if grid.xur < east or grid.yur < north:
        sys.exit("error: computed wind grid does not cover the landscape")
    return grid


def parse_start(text):
    """'MM DD HHMM' -> (month, day, hhmm)."""
    parts = re.split(r"[\s:/-]+", text.strip())
    if len(parts) != 3:
        sys.exit("error: --start must look like \"04 26 1200\" (month day HHMM), got %r" % text)
    try:
        mo, dy, hhmm = int(parts[0]), int(parts[1]), int(parts[2])
    except ValueError:
        sys.exit("error: --start must be numeric, got %r" % text)
    if not (1 <= mo <= 12 and 1 <= dy <= 31):
        sys.exit("error: --start month/day out of range: %r" % text)
    if not (0 <= hhmm <= 2359) or hhmm % 100 > 59:
        sys.exit("error: --start time must be HHMM in 0000-2359, got %r" % text)
    return mo, dy, hhmm


def read_input_file(path):
    """Pull FARSITE_START_TIME / FARSITE_END_TIME out of a FARSITE .input file."""
    start = end = None
    with open(path) as fh:
        for line in fh:
            m = re.match(r"\s*FARSITE_(START|END)_TIME:\s*(\d+)\s+(\d+)\s+(\d+)", line)
            if m:
                val = (int(m.group(2)), int(m.group(3)), int(m.group(4)))
                if m.group(1) == "START":
                    start = val
                else:
                    end = val
    if not start:
        sys.exit("error: no FARSITE_START_TIME in %s" % path)
    return start, end


# February is always 28 here, matching FARSITE: ConvertActualTimeToSimtime
# builds its reference from GetJulianDays(month) + day and adds a flat 365 on a
# year wrap, so the model itself has no leap-year concept. Keeping the same
# assumption avoids a one-day offset between the .atm records and sim time.
DAYS_IN_MONTH = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]


def step_times(month, day, hhmm, count, interval_min):
    """Yield (month, day, hhmm) every interval_min, rolling over days."""
    minutes = (hhmm // 100) * 60 + (hhmm % 100)
    for _ in range(count):
        yield month, day, (minutes // 60) * 100 + (minutes % 60)
        minutes += interval_min
        while minutes >= 1440:
            minutes -= 1440
            day += 1
            if day > DAYS_IN_MONTH[month - 1]:
                day = 1
                month += 1
                if month > 12:
                    month = 1


def minutes_between(a, b):
    """Rough elapsed minutes between two (month, day, hhmm) tuples."""
    def abs_min(t):
        mo, dy, hhmm = t
        doy = sum(DAYS_IN_MONTH[: mo - 1]) + dy
        return doy * 1440 + (hhmm // 100) * 60 + (hhmm % 100)
    return abs_min(b) - abs_min(a)


def smooth_field(grid, scale_cells, rng, sd):
    """Spatially correlated Gaussian field via a coarse grid + bilinear upsample.

    scale_cells is the correlation length. Values are drawn on a coarse lattice
    of that spacing and interpolated, so neighbouring cells stay similar --
    closer to a real wind field than per-cell white noise.

    Bilinear blending of independent draws shrinks the variance: averaged over a
    cell the four weights satisfy E[sum w^2] = (2/3)^2 = 4/9, so the realised
    standard deviation would come out at 2/3 of what was asked for. The coarse
    draws are inflated by the reciprocal so that the field actually delivers the
    requested sd.
    """
    if sd == 0.0:
        return None
    step = max(1, int(scale_cells))
    ccols = grid.ncols // step + 2
    crows = grid.nrows // step + 2
    coarse_sd = sd * 1.5 if step > 1 else sd
    coarse = [[rng.gauss(0.0, coarse_sd) for _ in range(ccols)] for _ in range(crows)]

    def at(col, row):
        fx, fy = col / float(step), row / float(step)
        x0, y0 = int(fx), int(fy)
        tx, ty = fx - x0, fy - y0
        c00, c10 = coarse[y0][x0], coarse[y0][x0 + 1]
        c01, c11 = coarse[y0 + 1][x0], coarse[y0 + 1][x0 + 1]
        return ((c00 * (1 - tx) + c10 * tx) * (1 - ty)
                + (c01 * (1 - tx) + c11 * tx) * ty)
    return at


def write_grid(path, grid, value_fn, fmt="%.2f"):
    """Write an ESRI ASCII grid; value_fn(col, row) -> float. Row 0 is the north edge."""
    with open(path, "w") as fh:
        fh.write(grid.header())
        for row in range(grid.nrows):
            fh.write("\t".join(fmt % value_fn(c, row) for c in range(grid.ncols)))
            fh.write("\n")


def load_hrrr(grid, when, args):
    """STUB: build (speed, direction) fields from HRRR forecast output.

    Not implemented. Filling this in needs, roughly:

      1. Locate the HRRR file covering `when` (GRIB2, hourly; typically
         hrrr.t{cycle}z.wrfsfcf{fhr}.grib2) under args.hrrr_dir.
      2. Read the 10 m wind components UGRD:10 m above ground and
         VGRD:10 m above ground -- e.g. via GDAL's GRIB driver and band
         metadata, or cfgrib/pygrib.
      3. Rotate the components from HRRR's Lambert Conformal grid-relative
         frame to true north. HRRR distributes grid-relative winds; skipping
         this rotation biases direction by the local grid convergence, up to
         ~20 degrees at the edges of CONUS.
      4. Reproject/resample onto `grid` (which is in the landscape's CRS) --
         bilinear on u and v separately, never on speed/direction, because
         direction wraps at 360.
      5. Convert: speed = hypot(u, v); direction_from =
         (270 - degrees(atan2(v, u))) mod 360.
      6. Return two callables (col, row) -> value, in m/s and degrees-from,
         matching the signature the constant/noisy paths use.

    Note step 5 yields a 10 m wind, so the output should stay METRIC and let
    FARSITE apply its own 10 m -> 20 ft reduction.
    """
    raise NotImplementedError(
        "HRRR ingest is not implemented yet. load_hrrr() in this file documents "
        "what it needs; --mode constant and --mode noisy work today.")


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ref = p.add_mutually_exclusive_group(required=True)
    ref.add_argument("--lcp", help="landscape file to match (header read directly, no GDAL)")
    ref.add_argument("--like", help="any GDAL-readable raster to match instead")
    p.add_argument("-o", "--outdir", required=True, help="output directory")
    p.add_argument("--prefix", default="wind", help="output file prefix (default: wind)")
    p.add_argument("--mode", choices=["constant", "noisy", "hrrr"], default="constant")

    p.add_argument("--speed", type=float, default=5.0,
                   help="wind speed in --speed-units (default: 5)")
    p.add_argument("--speed-units", choices=["mps", "kph", "mph"], default="mps",
                   help="units of --speed and --noise-speed (default: mps)")
    p.add_argument("--direction", type=float, default=270.0,
                   help="azimuth the wind blows FROM, degrees (default: 270, "
                        "a westerly blowing west-to-east)")
    p.add_argument("--toward", type=float,
                   help="azimuth the wind blows TOWARD; alternative to --direction")

    p.add_argument("--start", help='first time step, "MM DD HHMM" (e.g. "04 26 1200")')
    p.add_argument("--input-file",
                   help="FARSITE .input file to take the start (and end) time from")
    p.add_argument("--hours", type=float, default=6.0,
                   help="duration to cover, hours (default: 6; ignored if "
                        "--input-file supplies an end time)")
    p.add_argument("--interval-minutes", type=int, default=60,
                   help="time step in minutes; sub-hourly is written as HHMM "
                        "(default: 60)")

    p.add_argument("--units", choices=["metric", "english"], default="metric",
                   help="units keyword written into the .atm. metric = km/h at "
                        "10 m (FARSITE reduces to 20 ft); english = mph at 20 ft "
                        "(default: metric)")
    p.add_argument("--resolution", type=float,
                   help="wind grid cell size (default: same as the reference)")
    p.add_argument("--pad-cells", type=int, default=1,
                   help="cells of padding beyond the landscape (default: 1)")

    p.add_argument("--noise-speed", type=float, default=1.0,
                   help="speed noise standard deviation, in --speed-units "
                        "(mode=noisy; default: 1.0)")
    p.add_argument("--noise-direction", type=float, default=15.0,
                   help="direction noise standard deviation, degrees "
                        "(mode=noisy; default: 15)")
    p.add_argument("--noise-scale", type=float, default=10.0,
                   help="spatial correlation length in cells; 0 gives per-cell "
                        "white noise (mode=noisy; default: 10)")
    p.add_argument("--noise-seed", type=int, default=0,
                   help="RNG seed for reproducible noise (default: 0)")
    p.add_argument("--hrrr-dir", help="directory of HRRR GRIB2 files (mode=hrrr)")
    args = p.parse_args()

    if args.mode == "hrrr":
        try:
            load_hrrr(None, None, args)
        except NotImplementedError as e:
            sys.exit("error: %s" % e)

    # --- time window ---
    end_time = None
    if args.input_file:
        start_t, end_t = read_input_file(args.input_file)
        month, day, hhmm = start_t
        end_time = end_t
        print("time window from %s: start %02d %02d %04d%s"
              % (args.input_file, month, day, hhmm,
                 (", end %02d %02d %04d" % end_t) if end_t else ""))
    elif args.start:
        month, day, hhmm = parse_start(args.start)
    else:
        p.error("give --start or --input-file")

    if end_time:
        span = minutes_between((month, day, hhmm), end_time)
        if span <= 0:
            sys.exit("error: end time is not after start time")
        # One extra step so the last grid covers the end of the simulation.
        nsteps = int(math.ceil(span / float(args.interval_minutes))) + 1
    else:
        nsteps = int(math.ceil(args.hours * 60.0 / args.interval_minutes)) + 1
    if args.interval_minutes <= 0:
        sys.exit("error: --interval-minutes must be positive")

    # --- direction ---
    if args.toward is not None:
        direction = (args.toward + 180.0) % 360.0
        print("direction: blowing toward %g deg -> FROM %g deg"
              % (args.toward % 360.0, direction))
    else:
        direction = args.direction % 360.0

    # --- speed conversion ---
    to_mps = {"mps": 1.0, "kph": KPH_TO_MPS, "mph": MPH_TO_MPS}[args.speed_units]
    speed_mps = args.speed * to_mps
    noise_mps = args.noise_speed * to_mps
    if speed_mps < 0:
        sys.exit("error: --speed must not be negative")

    if args.units == "metric":
        out_speed = speed_mps * MPS_TO_KPH
        out_noise = noise_mps * MPS_TO_KPH
        unit_note = "km/h at 10 m"
        effective_mph = out_speed * TEN_M_TO_TWENTY_FT * (MPS_TO_MPH / MPS_TO_KPH)
    else:
        out_speed = speed_mps * MPS_TO_MPH
        out_noise = noise_mps * MPS_TO_MPH
        unit_note = "mph at 20 ft"
        effective_mph = out_speed

    grid = build_grid(args)
    print("wind grid: %s" % grid)
    print("speed:     %g %s -> %.4f %s" % (args.speed, args.speed_units, out_speed, unit_note))
    print("           FARSITE will use ~%.2f mph at 20 ft" % effective_mph)
    print("direction: %g deg FROM (blowing toward %g deg)" % (direction, (direction + 180) % 360))
    print("steps:     %d every %d min%s" % (nsteps, args.interval_minutes,
          "" if args.interval_minutes == 60 else " (sub-hourly)"))

    os.makedirs(args.outdir, exist_ok=True)
    times = list(step_times(month, day, hhmm, nsteps, args.interval_minutes))

    atm_lines = []
    if args.mode == "constant":
        # The field never changes, so write one pair and reference it from every
        # .atm record rather than duplicating identical grids.
        vel = os.path.join(args.outdir, "%s_vel.asc" % args.prefix)
        ang = os.path.join(args.outdir, "%s_ang.asc" % args.prefix)
        write_grid(vel, grid, lambda c, r: out_speed, "%.2f")
        write_grid(ang, grid, lambda c, r: direction, "%.0f")
        print("\nwrote %s" % os.path.basename(vel))
        print("wrote %s" % os.path.basename(ang))
        for mo, dy, hm in times:
            atm_lines.append((mo, dy, hm, vel, ang))
    else:
        rng = random.Random(args.noise_seed)
        for i, (mo, dy, hm) in enumerate(times):
            vel = os.path.join(args.outdir, "%s_%02d%02d_%04d_vel.asc" % (args.prefix, mo, dy, hm))
            ang = os.path.join(args.outdir, "%s_%02d%02d_%04d_ang.asc" % (args.prefix, mo, dy, hm))
            if args.noise_scale and args.noise_scale > 0:
                sfield = smooth_field(grid, args.noise_scale, rng, out_noise)
                dfield = smooth_field(grid, args.noise_scale, rng, args.noise_direction)
                spd_fn = (lambda c, r, f=sfield: max(0.0, out_speed + f(c, r))) if sfield \
                    else (lambda c, r: out_speed)
                dir_fn = (lambda c, r, f=dfield: (direction + f(c, r)) % 360.0) if dfield \
                    else (lambda c, r: direction)
            else:
                spd_fn = lambda c, r: max(0.0, rng.gauss(out_speed, out_noise))
                dir_fn = lambda c, r: (rng.gauss(direction, args.noise_direction)) % 360.0
            write_grid(vel, grid, spd_fn, "%.2f")
            write_grid(ang, grid, dir_fn, "%.0f")
            atm_lines.append((mo, dy, hm, vel, ang))
        print("\nwrote %d grid pairs" % len(times))

    atm_path = os.path.join(args.outdir, "%s.atm" % args.prefix)
    with open(atm_path, "w") as fh:
        # The keyword must be bare on its own line and upper case: the reader
        # compares it with strcmp, so "Metric" would be ignored and the grids
        # silently read as mph at 20 ft.
        fh.write("%s\n" % args.units.upper())
        for mo, dy, hm, vel, ang in atm_lines:
            # Absolute paths: FARSITE resolves grid paths against its own working
            # directory, not the .atm's location (CWindGrids::Create only chdir's
            # for Windows-style backslash paths), so a relative path here breaks
            # as soon as the binary is run from anywhere else.
            fh.write("%02d %02d %04d %s %s\n"
                     % (mo, dy, hm, os.path.abspath(vel), os.path.abspath(ang)))
    print("wrote %s  (%d records)" % (os.path.basename(atm_path), len(atm_lines)))

    print("\nAdd this to your FARSITE .input file:")
    print("  FARSITE_ATM_FILE: %s" % os.path.abspath(atm_path))
    print("The first record is %02d %02d %04d; FARSITE requires it at or before "
          "FARSITE_START_TIME." % times[0])
    return 0


if __name__ == "__main__":
    sys.exit(main())
