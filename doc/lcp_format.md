# The FARSITE landscape file (.lcp) format

Written for: anyone generating `.lcp` files for this version of FARSITE, or
maintaining its landscape reader.

Reverse-engineered from `src/Farsite5.cpp` (`ReadHeader`, `CellData`),
`src/Farsite5.h` (`headdata`), `src/fsxwenvt2.cpp` (unit conversion) and
`src/gridthem.cpp` (statistics), and checked against real files in `examples/`.
LANDFIRE stopped distributing `.lcp` in 2024, so anything you feed FARSITE now
has to be built by a converter -- see `scripts/tif2lcp.py`.

## Layout

A fixed **7316-byte** header, then the raster.

> `sizeof(headdata)` in memory is **7328**, not 7316, because of alignment
> padding after `latitude`. The on-disk layout has no padding. This is why
> `ReadHeader` reads field by field and why the single-`fread` attempt commented
> out at `Farsite5.cpp:717` would have been wrong. A writer must serialise field
> by field or pack the struct.

All integers are little-endian; `int32` is 4 bytes, `int16` 2, `double` 8.

### Header

| Offset | Type | Field | Notes |
|---|---|---|---|
| 0 | int32 | `CrownFuels` | 20 = absent, 21 = present |
| 4 | int32 | `GroundFuels` | 20 = absent, 21 = present |
| 8 | int32 | `latitude` | whole degrees; drives solar radiation |
| 12 | double×4 | `loeast, hieast, lonorth, hinorth` | **ignored on read** -- `ReadHeader:800-803` overwrites them from the UTM corners |
| 44 | — | 10 theme stat blocks, 412 bytes each | order below; ends at 4164 |
| 4164 | int32×2 | `numeast, numnorth` | columns, rows |
| 4172 | double×4 | `EastUtm, WestUtm, NorthUtm, SouthUtm` | the authoritative extent |
| 4204 | int32 | `GridUnits` | 0 = metric, 1 = English |
| 4208 | double×2 | `XResol, YResol` | **advisory**: `ReadHeader:823` recomputes resolution as `(EastUtm-WestUtm)/numeast`. Keep them consistent |
| 4224 | int16×10 | `EUnits, SUnits, AUnits, FOptions, CUnits, HUnits, BUnits, PUnits, DUnits, WOptions` | see [Unit codes](#unit-codes) |
| 4244 | char[256]×10 | per-theme source filenames | informational |
| 6804 | char[512] | `Description` | informational |

Each theme stat block is `lo` (int32), `hi` (int32), `num` (int32),
`cats[100]` (int32). Theme order throughout:

`0 elevation, 1 slope, 2 aspect, 3 fuel model, 4 canopy cover,
5 canopy height, 6 canopy base height, 7 canopy bulk density,
8 duff, 9 coarse woody debris`

### Raster

Row-major from the **north** edge down, west to east, interleaved: each cell is
`NumVals` consecutive `int16`s in theme order.

| `CrownFuels` | `GroundFuels` | `NumVals` | Themes stored |
|---|---|---|---|
| 20 | 20 | 5 | basic only |
| 20 | 21 | 7 | basic + duff/woody |
| 21 | 20 | 8 | basic + crown |
| 21 | 21 | 10 | all |

Total file size = `7316 + numeast * numnorth * NumVals * 2`.

## The statistics block is not decoration

`num` is the count of distinct values, or **-1** when there are more than 99.
The category list is **1-based**: slot 0 is unused, values occupy slots
`1..num` ascending. `lo`/`hi` exclude nodata.

Confirmed by the consumer loop at `Farsite5.cpp:9729`
(`for (i = 1; i <= grid->NumAllCats[3]; i++)`) and by every real file in
`examples/`. `LandscapeTheme::SortCats` (`gridthem.cpp:911`) emits exactly this
layout, so it is the de facto reference implementation.

**Getting this wrong fails silently and catastrophically.** `Far_Cond.cpp:125`
feeds `fuels[]`/`numfuel` into the fuel-moisture conditioning code, which keeps
only values that `GetFuelConversion` maps into `(0, 257)`. If the list is
malformed, no moisture keys get built, `NumFuels` ends at 0, and
`FE2::GetMx` (`FMC_FE22.cpp:48`) returns **0.0 for every 10-hr, 100-hr and
1000-hr moisture** -- the whole landscape burns bone dry, with nothing logged.

Not every field is read, though. Measured by which ones change results:

| Fields | Read by | Matters? |
|---|---|---|
| `fuels[]`, `numfuel` | `Set_FuelModel` → `AllocFuels` | **yes -- catastrophic** |
| `loelev`, `hielev`, `EUnits` | `AllocElevations`, `GetMx` | **yes** (band count and quantisation phase) |
| `loslope`, `hislope`, `SUnits` | `AllocSlopes` (clamps to 50°) | yes |
| `locover`, `hicover`, `CUnits` | `AllocCovers` | yes |
| `woodies[]`, `numwoody` | `Set_Woody` | no (`Far_Cond.cpp:118` says unused) |
| aspect/height/base/density/duff min-max, `lofuel`/`hifuel` | nothing | no |

If a file's statistics block is suspect, run FARSITE with
`--recompute-lcp-stats`, which rebuilds the whole block from the raster via
`LandscapeTheme::AnalyzeStats()` instead of trusting the header.

## Unit codes

Codes are recorded, **not applied** -- a wrong code silently reinterprets every
value. Semantics from `fsxwenvt2.cpp`.

| Field | Values | LANDFIRE |
|---|---|---|
| `EUnits` | 0 = metres, 1 = feet | 0 |
| `SUnits` | 0 = degrees, 1 = percent | 0 |
| `AUnits` | 0 = GRASS categories 1-25, 1 = GRASS degrees CCW-from-east, 2 = ArcInfo azimuth | **2** |
| `FOptions` | 0 = none, 1 = custom models, 2 = conversion file, 3 = both | 0 |
| `CUnits` | 0 = classes 1-4/99, 1 = percent | **1** |
| `HUnits`, `BUnits` | 1 = m, 2 = ft, 3 = m×10, 4 = ft×10 | 3 |
| `PUnits` | 1 = kg/m³, 2 = lb/ft³, 3 = kg/m³×100, 4 = lb/ft³×1000 | 3 |
| `DUnits` | 0 = class, 1 = tons/acre×10, 2 = Mg/ha×10 | — |
| `WOptions` | 0 = class code | — |

### Scaling and nodata traps

- **Canopy height, canopy base height: divided by 10 unconditionally.** The unit
  code only selects feet vs metres. Pixels must already be ×10.
- **Bulk density: divided by 100 unconditionally** → kg/m³.
- **Duff: divided by 10**; `DUnits == 2` then converts tons/acre → Mg/ha.
- `AUnits == 0` with the value 25 forces `slope = 0` as well as aspect.
- `CUnits == 0` maps classes to 10/30/60/75 %, `99` → 0 %, and **anything else
  to 0 %** via the `default:` branch. Feeding percent data with `CUnits = 0`
  therefore zeroes all canopy cover.
- **Silent rescue heuristics.** `HeightConvert`/`BaseConvert` contain
  `if (val > 100.0) val /= 10.0; // probably got wrong units`, and
  `DensityConvert` has `if (density > 1.0) density /= 100.0`. These fire on
  *legitimate* data (a 120 m canopy, CBD above 1.0 kg/m³) and will partially
  mask a 10× scaling error, so a mis-scaled file can produce
  plausible-but-wrong output. Test with values straddling the thresholds.
- **Nodata is `-9999` by convention and handled inconsistently per theme**:
  elevation, slope and aspect compare `== -9999` exactly, while canopy height,
  base, density and duff test `>= 0` and substitute defaults. Fuel goes through
  `GetFuelConversion`, which maps anything outside `[0, 257)` to -1
  (unburnable). `-9999` belongs in the category list but **must not** appear in
  `lo<theme>`.

## Geometry constraints

One grid for all themes: north-up, square, uniformly spaced cells.
`GetCellPosition` (`Farsite5.cpp:906`) assumes it, and there is **no CRS
handling anywhere** in this codebase -- `ConvertUtmToEastingOffset` and friends
(`Farsite5.cpp:958-1000`) are identity functions, and output grids are written
without a `.prj`. Coordinates are whatever the input used; keep every input and
your ignition/barrier shapefiles in one projection.

When resampling to a common grid, use **nearest neighbour** for fuel model,
categorical aspect and class-coded cover: interpolating class codes invents
values that do not exist, which then silently become unburnable.
