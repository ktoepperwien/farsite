# Building a FARSITE case from scratch

Written for: anyone setting up a new fire-spread run in this fork.

FARSITE needs four things: a landscape file, an ignition, weather, and an input
file tying them together. The scripts in `scripts/` generate each one and check
the constraints that otherwise fail silently. `examples/Dish/run_dish.sh` runs
the whole chain.

```
LANDFIRE GeoTIFFs ──crop_landfire.py──> per-theme .tif ──tif2lcp.py──> .lcp
                                                                        │
                          ┌─────────────────────────────────────────────┤
                          │                     │                       │
                   make_wind.py          make_ignition.py         make_input.py
                          │                     │                       │
                    wind.atm + grids       ignition.shp          case.input + .raws
                          └─────────────────────┴───────────────────────┘
                                              │
                                      TestFARSITE <runfile>
```

## Quick start

> **Note:** the `examples/Dish/` case used throughout this document is **not
> tracked in git** — it is ~36 MB of mostly regenerable artefacts, and its outputs
> are useless as a cross-platform reference because these floating-point-sensitive
> cases diverge by architecture. Recreate it with `crop_landfire.py` from LANDFIRE
> CONUS mosaics (section 1), or substitute any landscape of your own. The tracked
> `examples/Panther`, `examples/cougarCreek` and `examples/flatland` cases work out
> of the box and are what `tests/run_regression.sh` uses.

```
./examples/Dish/run_dish.sh
```

Everything is overridable by environment variable:

```
WIND_MODE=noisy WIND_SPEED=8 WIND_DIR=225 HOURS=12 SCENARIO=very-dry \
  ./examples/Dish/run_dish.sh
```

---

## 1. Landscape

See `doc/lcp_format.md` for the format itself. In short:

```
# crop LANDFIRE CONUS mosaics to a window around a fire centroid
python3 scripts/crop_landfire.py \
    --available-manifest available.csv --training-manifest training.csv \
    --fire-id 0 --size-cells 512 -o examples/Dish/

# stack the cropped GeoTIFFs into a landscape file
examples/Dish/0_Dish/make_lcp.sh
```

`--size-cells` is in multiples of the native 30 m cell, so 512 gives a 15.36 km
box. `--resolution` resamples; the window is snapped to the source grid, so at
native resolution the crop is an exact pixel subset.

Themes LANDFIRE does not ship (canopy cover, height, base height, bulk density)
are written as constants — zero by default, meaning no canopy. Duff and coarse
woody debris are omitted entirely rather than faked, so the landscape declares
`GroundFuels = absent`.

**Do not use `gdal_translate -of LCP`.** Its pixel data is right but its header
statistics block is not, and a wrong fuel-model category list makes every
conditioned dead fuel moisture come back `0.0` with nothing logged. If you have
a landscape file from elsewhere and results look wrong, run FARSITE with
`--recompute-lcp-stats` to rebuild that block from the raster.

## 2. Wind

```
python3 scripts/make_wind.py --lcp case.lcp -o case/wind \
    --start "04 26 1200" --hours 6 --speed 5 --direction 270
```

| Option | Default | Notes |
|---|---|---|
| `--mode` | `constant` | `constant`, `noisy`, or `hrrr` (stub) |
| `--speed` | `5` | in `--speed-units` (`mps`, `kph`, `mph`) |
| `--direction` | `270` | azimuth the wind blows **from**; 270 is a westerly blowing west-to-east |
| `--toward` | — | name the direction the wind is heading instead |
| `--hours` / `--input-file` | `6` | duration, or take the window from a `.input` |
| `--interval-minutes` | `60` | sub-hourly is supported and written as `HHMM` |
| `--units` | `metric` | see below |
| `--noise-speed` / `--noise-direction` | `1.0` / `15` | std dev; only with `--mode noisy` |
| `--noise-scale` | `10` | spatial correlation length in cells; `0` gives white noise |
| `--resolution` | LCP cell | wind grid need not match the landscape |
| `--pad-cells` | `1` | padding beyond the landscape |

### Wind units and direction, which are easy to get wrong

The `.atm` units keyword decides how FARSITE reads the speed grids:

- `METRIC` — **km/h at 10 m**. FARSITE divides by 1.15 for a 20 ft wind, then
  converts to mph.
- `ENGLISH` — **mph at 20 ft**, used as-is. This is the default when the keyword
  is absent.

Neither is m/s. `--speed 5` (m/s) is written as 18.0 km/h under the default
metric output, and FARSITE ends up using ~9.73 mph at 20 ft. Metric is the
default so FARSITE applies its own height reduction rather than the script
guessing at it. The keyword must be bare and upper case — the reader uses
`strcmp`, so `Metric` is ignored and the grids get read as mph.

Direction is the meteorological "from" convention. This was confirmed by
experiment, not by reading the trig: wind from 270° drove the fire to azimuth
94°, and `--toward 0` drove it to azimuth 8°.

Two constraints FARSITE enforces: the **first `.atm` record must be at or
before `FARSITE_START_TIME`**, and the **wind grid must cover the landscape**.
`--input-file` handles the first automatically. Grid paths are written absolute,
because FARSITE resolves them against its own working directory rather than the
`.atm`'s location.

## 3. Ignition

```
python3 scripts/make_ignition.py --lcp case.lcp -o case/ignition.shp --radius 200
```

| Option | Default | Notes |
|---|---|---|
| `--mode` | `circle` | `circle`, or `check` to validate an existing shapefile |
| `--center X Y` | landscape centre | in landscape coordinates |
| `--radius` | `200` | metres |
| `--vertices` | `32` | |
| `--winding` | `ccw` | `ccw` matches FARSITE natively; either works |
| `--input` | — | shapefile to validate (`--mode check`) |
| `--point-to-circle` | off | convert point ignitions to `--radius` circles |
| `--min-edge-cells` | `2` | required clearance from the landscape edge |
| `--force` | off | write despite validation problems |

Both modes report edge clearance, the cells covered, and the fuel under the
ignition, then refuse to write on a real problem. The two silent failure modes
are landing outside the landscape — everything off-grid reads as fuel `-9999`,
which converts to unburnable — and landing entirely on non-burnable fuel.

Geometry types, per `IgnitionFile::ShapeInput`:

- **Polygon** — an area ignition. Winding is normalised internally (`arp()`
  computes a signed area and reverses a negative ring), verified: clockwise and
  counter-clockwise inputs give byte-identical output.
- **PolyLine** — a line source. Needs the `fsxwignt.cpp:191` fix; on an unfixed
  build the simulation never terminates. See `tests/test_line_ignition.sh`.
- **Point** — FARSITE expands each point to a 10-vertex circle of radius
  `startsize`, hard-coded to **1 map unit**. On a 30 m landscape that is far
  below one cell, so prefer a polygon or use `--point-to-circle`.

There is **no reprojection anywhere in FARSITE**. The ignition, the wind grids
and the landscape must already share one coordinate system.

## 4. Input file, fuel moisture and weather

```
python3 scripts/make_input.py --lcp case.lcp \
    --atm case/wind/wind.atm --make-raws case/weather.raws \
    -o case/case.input --start "04 26 1200" --hours 6 \
    --moisture-mode per-model --moisture-noise 1.5
```

### Is fuel moisture a scalar or an array?

**An array of per-fuel-model records, each holding five integer scalars.**

```
FUEL_MOISTURES_DATA: <n>
<fuel_model> <1hr> <10hr> <100hr> <live_herb> <live_woody>
```

Three dead classes by time lag (1, 10, 100 hour) plus two live classes
(herbaceous, woody), all whole percent. There is no 1000-hour column here even
though the moisture model tracks one internally.

Fuel model `0` is the default record. `SetAllMoistures` (Farsite5.cpp:1930)
writes it into all 257 model slots, so a single model-0 line satisfies every
fuel model on the landscape; any later record overrides that one model. Values
are clamped to a floor of 2 and validated against 2–300 (`e_GFMlow` /
`e_GLFMup`).

It is **not** a spatial field. You cannot set moisture per cell. Spatial and
temporal variation comes from the conditioning pass, which builds moisture over
elevation × slope × aspect × canopy-cover bands driven by the weather stream,
using these numbers as the starting point.

So per-class random noise only means anything *across fuel models*, which is
what `--moisture-mode per-model` does: one record per fuel model actually
present on the landscape, each perturbed. On the Dish landscape that is 22
burnable models out of 26 present. Live classes get 4× the dead-class noise,
since they vary over a much wider range.

| Option | Default | Notes |
|---|---|---|
| `--scenario` | `dry` | `very-dry` 3/4/5/30/60, `dry` 4/5/7/60/90, `moderate` 6/8/10/90/120, `wet` 10/12/14/120/150 |
| `--moisture H1 H10 H100 HERB WOODY` | — | explicit, overrides `--scenario` |
| `--moisture-mode` | `default` | `default` (one model-0 record) or `per-model` |
| `--moisture-noise` | `0` | std dev in percentage points on the dead classes |
| `--timestep` | `30` | minutes |
| `--distance-res` / `--perimeter-res` | cell / 2×cell | |
| `--spot-probability` | `0` | 0 disables spotting |
| `--crown-fire-method` | `Finney` | or `ScottReinhardt` |

### Weather is also required

Fuel conditioning needs temperature, humidity and precipitation, so a weather
stream is mandatory even when gridded wind supplies the wind. `--make-raws`
writes a uniform RAWS stream with `--raws-lead-hours` (default 24) of lead-in
before the start so the fuels have something to condition against.

```
RAWS_ELEVATION: 500
RAWS_UNITS: English
RAWS: 31
2013 4 25 1200 77 25 0.00 3 270 25
#year month day HHMM temp(F) humid(%) precip(in) wspeed(mph) wdir(from) cloud(%)
```

`RAWS_UNITS` is compared with a case-sensitive `strcmp` on this branch, so it
must be exactly `English` or `Metric` — `ENGLISH` fails validation. When an
`.atm` file is present the RAWS wind columns are ignored for spread, but they
still have to be there.

The repo's own history records that gridded winds worked with `.atm` **plus** a
RAWS file but not with `.atm` alone, because RAWS triggers the interpolation
that extends weather across the run.

## 5. Running

The command file takes one run per line, six tokens:

```
<landscape.lcp> <case.input> <ignition.shp> <barrier.shp|0> <output_base> <type>
```

`type`: `0` both ASCII and binary, `1` ASCII grids, `2` FlamMap binary, `4`
ASCII plus hourly fuel-moisture maps.

```
./src/TestFARSITE examples/Dish/0_Dish/run_dish.txt
```

Optional flag: `--recompute-lcp-stats` rebuilds the landscape header's
statistics from the raster instead of trusting it. Off by default so
well-formed files reproduce exactly.

## 6. Plotting the result

```
/Users/karl/anaconda3/envs/geo/bin/python scripts/plot_arrival.py \
    --case examples/Dish/0_Dish/output/dish -o dish_arrival.png \
    --contour-labels --perimeters examples/Dish/0_Dish/output/dish_Perimeters.shp
```

Draws the arrival-time grid with `pcolormesh` in EPSG:5070 over ESRI World
Imagery (contextily + xyzservices), with gray contour lines of the arrival
front. This needs numpy, matplotlib, contextily and xyzservices, which the
repo's other scripts do not — on this machine the `geo` conda environment has
them, hence the explicit interpreter above.

| Option | Default | Notes |
|---|---|---|
| `--case` / `--grid` | — | output base (appends `_ArrivalTime.asc`), or a grid directly |
| `--crs` | `EPSG:5070` | CRS of the grid, used for the basemap fetch |
| `--units` | `hours` | FARSITE writes minutes; `hours` converts |
| `--cmap` | `inferno` | sequential ramp; rainbow ramps are refused |
| `--cmap-floor` | `0.18` | skip this much of the ramp's dark end |
| `--alpha` | `0.80` | fill opacity over the imagery |
| `--contour-levels` | `6` | `0` disables |
| `--contour-labels` | off | inline labels, thinned by `--max-labelled-levels` |
| `--perimeters` | — | overlay a perimeter shapefile |
| `--no-basemap` | off | skip the tile fetch (works offline) |

Three choices worth knowing about:

**The fill is a heat ramp on purpose.** Arrival time is a magnitude, so it takes
a sequential ramp. Multi-hue sequential is normally wrong, but "semantic heat"
is a sanctioned exception and it ships with a colorbar as that exception
requires. Rainbow ramps (`jet`, `turbo`, `rainbow`, `hsv`) are rejected outright
— they are not perceptually uniform and invent structure that is not in the data.

**The contours are gray on purpose.** The fill already owns colour for
magnitude; the contours are a second, redundant reading of the same field, so
they stay neutral rather than competing for the same channel. Their gray ramps
light-to-dark with time so the front direction reads without a second legend.

**Zero is a burned value.** FARSITE writes arrival time `0` for the cells of the
initial ignition perimeter, so the mask excludes `< 0` rather than `<= 0`.
Masking zero punches a hole in the map exactly where the fire started. Cells
that remain unmasked-but-blank inside the burn are genuinely non-burnable fuel.

The default ramp's dark end is clipped (`--cmap-floor`) because heat ramps bottom
out at near-black, which disappears against dark satellite imagery.

## Troubleshooting

| Symptom | Likely cause |
|---|---|
| `Error opening output file` while loading inputs | a referenced file cannot be opened — usually a relative path in the `.atm`, or a missing RAWS |
| Nothing burns | ignition on non-burnable fuel, a too-wet scenario, or the ignition outside the landscape. `make_ignition.py` reports the fuel underneath |
| Run never terminates | polyline ignition on a build without the `fsxwignt.cpp:191` fix |
| Every dead fuel moisture is 0 | corrupt landscape statistics block; re-run with `--recompute-lcp-stats` |
| Wind seems ~1.6× too strong or weak | `.atm` units keyword mismatch: km/h-at-10 m versus mph-at-20 ft |
| Fire spreads the wrong way | direction is the azimuth the wind blows **from**, not toward |
| Wind grid rejected | first record after `FARSITE_START_TIME`, or grid smaller than the landscape |

## Reproducibility

Results are exactly reproducible on the platform that produced the committed
reference outputs, and only approximately reproducible elsewhere.

- **x86_64 Linux / gcc:** every case reproduces bit-for-bit — measured 207 of 207
  output files identical with gcc 12.4.0, at both `-O0` and `-O2`.
  `tests/run_regression.sh` therefore gates the whole suite strictly there by
  default.
- **Other platforms:** FARSITE's spread front is threshold-driven, so last-bit
  floating-point differences amplify into macroscopic ones. On Apple silicon,
  Panther's `cust` case still matches exactly but `test1177973`, `cougarCreek` and
  `flatland` differ by a few percent of burned area — and differ that much between
  `-O0` and `-O2` on the same machine with identical sources. The suite reports
  those without failing; `STRICT_ALL=1` forces it to gate on them anyway.

So a difference on x86_64 Linux is unambiguously a bug, while a difference on
another architecture usually is not. Check there first.

To compare against older unoptimised numbers, rebuild with
`make clean && make CXXFLAGS="-std=c++11 -g -Wall -DUNIX -Wno-deprecated"`.

Seeds worth pinning: `--noise-seed` (wind), `--moisture-seed`, and
`SPOTTING_SEED` in the input file.
