#!/usr/bin/env python3
"""Generate a FARSITE .input file (and optionally a RAWS weather stream).

Fuel moisture in FARSITE
------------------------
It is not a single scalar, and it is not a spatial grid either. The input is an
**array of per-fuel-model records, each holding five integer scalars**:

    FUEL_MOISTURES_DATA: <n>
    <fuel_model> <1hr> <10hr> <100hr> <live_herb> <live_woody>

Three dead classes by time-lag (1, 10, 100 hour), plus two live classes
(herbaceous and woody). All are whole percent. There is no 1000-hour class in
this input even though the moisture model tracks one internally.

Fuel model `0` is the default record (`e_DefFulMod`): `SetAllMoistures`
(Farsite5.cpp:1930) writes it into all 257 model slots, so one model-0 line
satisfies every fuel model on the landscape. Any further record overrides that
one model. Values are clamped to a floor of 2 and validated against 2..300
(`e_GFMlow`/`e_GLFMup` in icf_def.h).

So per-class random noise is only meaningful *across fuel models* -- which is
what --moisture-mode per-model does here. Spatial and temporal variation is not
yours to set: FARSITE derives it in the conditioning pass, building moisture
over elevation x slope x aspect x canopy-cover bands driven by the weather
stream, using these values as the starting point.

Example
-------
  python3 scripts/make_input.py \
      --lcp examples/Dish/0_Dish/0_Dish.lcp \
      --atm examples/Dish/0_Dish/wind/wind.atm \
      --make-raws examples/Dish/0_Dish/weather.raws \
      -o examples/Dish/0_Dish/dish.input \
      --start "04 26 1200" --hours 6 --moisture-mode per-model
"""

import argparse
import os
import random
import struct
import sys

LCP_HEADER_SIZE = 7316
LCP_OFF_CROWN = 0
LCP_OFF_FUELBLOCK = 44 + 3 * 412      # lofuel, hifuel, numfuel, fuels[100]
LCP_OFF_NUMEAST = 4164

# icf_def.h: e_GFMlow / e_GLFMup. Dead fuels are additionally capped at
# something physically sane -- the validator would accept 300% dead fuel
# moisture, which is nonsense.
MOIST_MIN, MOIST_MAX = 2, 300
DEAD_MAX, LIVE_MIN = 60, 30

# 1hr, 10hr, 100hr, live herbaceous, live woody -- whole percent.
SCENARIOS = {
    "very-dry": (3, 4, 5, 30, 60),
    "dry":      (4, 5, 7, 60, 90),
    "moderate": (6, 8, 10, 90, 120),
    "wet":      (10, 12, 14, 120, 150),
}

# Fuel models that cannot carry fire; no point emitting moisture for them.
NONBURNABLE = set(range(91, 100)) | {0, -9999}


def lcp_fuel_models(path):
    """Distinct fuel models on the landscape.

    Uses the header's category list when it is usable -- it is exactly this,
    already computed, 1-based with slot 0 unused. Falls back to scanning the
    raster when the theme has more than 99 distinct values (numfuel == -1).
    """
    with open(path, "rb") as fh:
        head = fh.read(LCP_HEADER_SIZE)
    if len(head) < LCP_HEADER_SIZE:
        sys.exit("error: %s is too short to be an .lcp" % path)
    lo, hi, num = struct.unpack_from("<3i", head, LCP_OFF_FUELBLOCK)
    cats = struct.unpack_from("<100i", head, LCP_OFF_FUELBLOCK + 12)
    if num > 0:
        return sorted({c for c in cats[1:num + 1]}), "header category list"

    crown, ground = struct.unpack_from("<2i", head, LCP_OFF_CROWN)
    numvals = {(20, 20): 5, (20, 21): 7, (21, 20): 8, (21, 21): 10}.get((crown, ground))
    if numvals is None:
        sys.exit("error: %s has CrownFuels=%d GroundFuels=%d" % (path, crown, ground))
    numeast, numnorth = struct.unpack_from("<2i", head, LCP_OFF_NUMEAST)
    found = set()
    with open(path, "rb") as fh:
        fh.seek(LCP_HEADER_SIZE)
        for _ in range(numnorth):
            raw = fh.read(numeast * numvals * 2)
            if len(raw) < numeast * numvals * 2:
                break
            row = struct.unpack("<%dh" % (numeast * numvals), raw)
            found.update(row[3::numvals])
    return sorted(found), "raster scan (header said >99 classes)"


def clamp(v, lo, hi):
    return max(lo, min(hi, int(round(v))))


def moisture_records(models, base, mode, noise, rng):
    """Build (model, 1h, 10h, 100h, herb, woody) records."""
    d1, d10, d100, lh, lw = base
    default = (0, clamp(d1, MOIST_MIN, DEAD_MAX), clamp(d10, MOIST_MIN, DEAD_MAX),
               clamp(d100, MOIST_MIN, DEAD_MAX), clamp(lh, LIVE_MIN, MOIST_MAX),
               clamp(lw, LIVE_MIN, MOIST_MAX))
    if mode == "default":
        return [default]

    recs = [default]        # keep model 0 so any model not listed is still covered
    for m in models:
        if m in NONBURNABLE:
            continue
        if noise > 0:
            vals = (clamp(rng.gauss(d1, noise), MOIST_MIN, DEAD_MAX),
                    clamp(rng.gauss(d10, noise), MOIST_MIN, DEAD_MAX),
                    clamp(rng.gauss(d100, noise), MOIST_MIN, DEAD_MAX),
                    clamp(rng.gauss(lh, noise * 4), LIVE_MIN, MOIST_MAX),
                    clamp(rng.gauss(lw, noise * 4), LIVE_MIN, MOIST_MAX))
        else:
            vals = (clamp(d1, MOIST_MIN, DEAD_MAX), clamp(d10, MOIST_MIN, DEAD_MAX),
                    clamp(d100, MOIST_MIN, DEAD_MAX), clamp(lh, LIVE_MIN, MOIST_MAX),
                    clamp(lw, LIVE_MIN, MOIST_MAX))
        recs.append((m,) + vals)
    return recs


def parse_start(text):
    parts = text.replace(":", " ").replace("/", " ").split()
    if len(parts) != 3:
        sys.exit('error: --start must look like "04 26 1200" (month day HHMM)')
    mo, dy, hhmm = (int(p) for p in parts)
    if not (1 <= mo <= 12 and 1 <= dy <= 31 and 0 <= hhmm <= 2359 and hhmm % 100 < 60):
        sys.exit("error: --start out of range: %r" % text)
    return mo, dy, hhmm


DAYS_IN_MONTH = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]


def add_hours(mo, dy, hhmm, hours):
    """Shift a (month, day, HHMM) by a signed number of hours.

    Negative shifts matter: the RAWS stream needs a conditioning lead-in that
    starts before the simulation, which routinely crosses back over midnight.
    February is 28 days here to match FARSITE's own leap-agnostic date maths.
    """
    minutes = (hhmm // 100) * 60 + hhmm % 100 + int(round(hours * 60))
    while minutes >= 1440:
        minutes -= 1440
        dy += 1
        if dy > DAYS_IN_MONTH[mo - 1]:
            dy = 1
            mo = 1 if mo == 12 else mo + 1
    while minutes < 0:
        minutes += 1440
        dy -= 1
        if dy < 1:
            mo = 12 if mo == 1 else mo - 1
            dy = DAYS_IN_MONTH[mo - 1]
    return mo, dy, (minutes // 60) * 100 + minutes % 60


def write_raws(path, start, hours, args):
    """Minimal RAWS stream covering the run, one record per hour.

    RAWS_UNITS is compared with a case-sensitive strcmp on this branch, so it
    has to be exactly "English" or "Metric" -- "ENGLISH" fails validation.
    English columns are: degF, %, inches, mph, degrees-from, % cloud.
    """
    mo, dy, hhmm = start
    # Cover the conditioning lead-in as well as the run itself: the moisture
    # model wants weather before ignition to condition the fuels.
    lead = args.raws_lead_hours
    cmo, cdy, chhmm = add_hours(mo, dy, hhmm, -lead)
    n = int(round(lead + hours)) + 1
    folder = os.path.dirname(os.path.abspath(path))
    if folder and not os.path.isdir(folder):
        os.makedirs(folder)
    with open(path, "w") as fh:
        fh.write("RAWS_ELEVATION: %d\n" % args.raws_elevation)
        fh.write("RAWS_UNITS: English\n")
        fh.write("RAWS: %d\n" % n)
        t = (cmo, cdy, chhmm)
        for _ in range(n):
            fh.write("%d %d %d %04d %d %d %.2f %d %d %d\n"
                     % (args.raws_year, t[0], t[1], t[2], args.temperature,
                        args.humidity, 0.0, args.raws_windspeed,
                        args.raws_winddir, args.cloud_cover))
            t = add_hours(t[0], t[1], t[2], 1)
    return n, (cmo, cdy, chhmm)


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--lcp", required=True, help="landscape file (read for its fuel models)")
    p.add_argument("-o", "--output", required=True, help="output .input path")
    p.add_argument("--start", default="04 26 1200", help='start "MM DD HHMM" (default: 04 26 1200)')
    p.add_argument("--hours", type=float, default=6.0, help="run duration, hours (default: 6)")
    p.add_argument("--timestep", type=int, default=30, help="FARSITE_TIMESTEP, minutes (default: 30)")
    p.add_argument("--distance-res", type=float, help="FARSITE_DISTANCE_RES (default: LCP cell size)")
    p.add_argument("--perimeter-res", type=float, help="FARSITE_PERIMETER_RES (default: 2x cell size)")

    p.add_argument("--scenario", choices=sorted(SCENARIOS), default="dry",
                   help="fuel moisture preset (default: dry)")
    p.add_argument("--moisture", nargs=5, type=int,
                   metavar=("H1", "H10", "H100", "HERB", "WOODY"),
                   help="explicit moisture percentages, overriding --scenario")
    p.add_argument("--moisture-mode", choices=["default", "per-model"], default="default",
                   help="'default' writes one model-0 record covering every fuel "
                        "model; 'per-model' also writes one record per fuel model "
                        "present on the landscape (default: default)")
    p.add_argument("--moisture-noise", type=float, default=0.0,
                   help="std dev of per-model noise on the dead classes, in "
                        "percentage points; live classes get 4x this. Only has an "
                        "effect with --moisture-mode per-model (default: 0)")
    p.add_argument("--moisture-seed", type=int, default=0, help="RNG seed (default: 0)")

    p.add_argument("--atm", help="gridded wind .atm file (FARSITE_ATM_FILE)")
    p.add_argument("--raws", help="existing RAWS weather stream to reference")
    p.add_argument("--make-raws", metavar="PATH",
                   help="generate a uniform RAWS stream at PATH and reference it")
    p.add_argument("--raws-lead-hours", type=float, default=24.0,
                   help="hours of weather before the start, for fuel conditioning "
                        "(default: 24)")
    p.add_argument("--raws-year", type=int, default=2013, help="year in the RAWS records")
    p.add_argument("--raws-elevation", type=int, default=500, help="RAWS_ELEVATION, feet")
    p.add_argument("--temperature", type=int, default=77, help="degF (default: 77)")
    p.add_argument("--humidity", type=int, default=25, help="%% (default: 25)")
    p.add_argument("--cloud-cover", type=int, default=25, help="%% (default: 25)")
    p.add_argument("--raws-windspeed", type=int, default=3,
                   help="mph in the RAWS stream. Ignored for spread when --atm "
                        "supplies gridded wind, but still required (default: 3)")
    p.add_argument("--raws-winddir", type=int, default=270,
                   help="degrees the wind blows FROM (default: 270)")

    p.add_argument("--foliar-moisture", type=int, default=100, help="%% (default: 100)")
    p.add_argument("--crown-fire-method", choices=["Finney", "ScottReinhardt"],
                   default="Finney")
    p.add_argument("--spot-probability", type=float, default=0.0,
                   help="0 disables spotting (default: 0)")
    p.add_argument("--spotting-seed", type=int, default=253)
    p.add_argument("--no-acceleration", action="store_true")
    args = p.parse_args()

    if not (args.raws or args.make_raws):
        sys.exit("error: FARSITE needs a weather stream for fuel conditioning. "
                 "Pass --raws <file> or --make-raws <path>.")
    if args.raws and args.make_raws:
        sys.exit("error: give either --raws or --make-raws, not both")

    # Cell size drives the sensible defaults for the resolution switches.
    with open(args.lcp, "rb") as fh:
        head = fh.read(LCP_HEADER_SIZE)
    numeast, numnorth = struct.unpack_from("<2i", head, LCP_OFF_NUMEAST)
    east, west, north, south = struct.unpack_from("<4d", head, 4172)
    cell = (east - west) / numeast
    dist_res = args.distance_res or cell
    perim_res = args.perimeter_res or 2 * cell

    models, how = lcp_fuel_models(args.lcp)
    burnable = [m for m in models if m not in NONBURNABLE]
    print("landscape: %s" % args.lcp)
    print("  %d x %d cells at %g" % (numeast, numnorth, cell))
    print("  fuel models (%s): %s" % (how, ", ".join(str(m) for m in models)))
    print("  burnable: %d of %d" % (len(burnable), len(models)))

    base = tuple(args.moisture) if args.moisture else SCENARIOS[args.scenario]
    rng = random.Random(args.moisture_seed)
    recs = moisture_records(burnable, base, args.moisture_mode, args.moisture_noise, rng)
    print("\nfuel moisture: %s, base %s (1h/10h/100h/herb/woody)"
          % (args.moisture_mode, "/".join(str(v) for v in base)))
    if args.moisture_mode == "per-model":
        print("  %d records (model 0 default + %d fuel models)%s"
              % (len(recs), len(recs) - 1,
                 ", noise sd %.1f" % args.moisture_noise if args.moisture_noise else ""))
    else:
        print("  1 record (model 0) -- SetAllMoistures applies it to every model")

    start = parse_start(args.start)
    end = add_hours(start[0], start[1], start[2], args.hours)
    print("\nrun: %02d %02d %04d -> %02d %02d %04d (%g h), timestep %d min"
          % (start + end + (args.hours, args.timestep)))

    raws_path = args.raws
    if args.make_raws:
        n, cstart = write_raws(args.make_raws, start, args.hours, args)
        raws_path = args.make_raws
        print("weather: wrote %s (%d hourly records from %02d %02d %04d, "
              "%g h of conditioning lead-in)"
              % (args.make_raws, n, cstart[0], cstart[1], cstart[2], args.raws_lead_hours))
    else:
        print("weather: %s" % raws_path)

    folder = os.path.dirname(os.path.abspath(args.output))
    if folder and not os.path.isdir(folder):
        os.makedirs(folder)
    L = []
    L.append("FARSITE INPUTS FILE VERSION 1.0")
    L.append("FARSITE_START_TIME: %02d %02d %04d" % start)
    L.append("FARSITE_END_TIME: %02d %02d %04d" % end)
    L.append("FARSITE_TIMESTEP: %d" % args.timestep)
    L.append("FARSITE_DISTANCE_RES: %.1f" % dist_res)
    L.append("FARSITE_PERIMETER_RES: %.1f" % perim_res)
    L.append("FARSITE_MIN_IGNITION_VERTEX_DISTANCE: %.1f" % (dist_res / 2.0))
    L.append("FARSITE_SPOT_GRID_RESOLUTION: %.1f" % (dist_res / 2.0))
    L.append("FARSITE_SPOT_PROBABILITY: %g" % args.spot_probability)
    L.append("FARSITE_SPOT_IGNITION_DELAY: 0")
    L.append("FARSITE_MINIMUM_SPOT_DISTANCE: %d" % int(round(cell)))
    L.append("FARSITE_ACCELERATION_ON: %d" % (0 if args.no_acceleration else 1))
    L.append("FARSITE_FILL_BARRIERS: 1")
    L.append("SPOTTING_SEED: %d" % args.spotting_seed)
    L.append("")
    L.append("# model 1hr 10hr 100hr live_herb live_woody   (whole percent)")
    L.append("# model 0 is the default and covers every fuel model.")
    L.append("FUEL_MOISTURES_DATA: %d" % len(recs))
    for r in recs:
        L.append("%d %d %d %d %d %d" % r)
    L.append("")
    if args.atm:
        L.append("# Gridded wind. Wind for spread comes from here, not from the")
        L.append("# RAWS wind columns; the weather stream is still needed for the")
        L.append("# temperature/humidity/precipitation that drives conditioning.")
        L.append("FARSITE_ATM_FILE: %s" % os.path.abspath(args.atm))
        L.append("")
    L.append("# year month day HHMM temp(F) humid(%) precip(in) wspeed(mph) wdir(deg from) cloud(%)")
    L.append("RAWS_FILE: %s" % os.path.abspath(raws_path))
    L.append("")
    L.append("FOLIAR_MOISTURE_CONTENT: %d" % args.foliar_moisture)
    L.append("CROWN_FIRE_METHOD: %s" % args.crown_fire_method)
    L.append("NUMBER_PROCESSORS: 1")
    L.append("")
    L.append("FLAMELENGTH:")
    L.append("SPREADRATE:")
    L.append("INTENSITY:")
    L.append("CROWNSTATE:")
    with open(args.output, "w") as fh:
        fh.write("\n".join(L) + "\n")

    print("\nwrote %s" % args.output)
    if args.atm:
        print("\nNote: FARSITE requires the first .atm record at or before")
        print("      FARSITE_START_TIME, and the wind grid must cover the landscape.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
