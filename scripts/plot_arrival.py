#!/usr/bin/env python3
"""Map a FARSITE fire arrival-time grid over satellite imagery.

Reads the ESRI ASCII grid FARSITE writes (`<base>_ArrivalTime.asc`), draws it
with pcolormesh in EPSG:5070, overlays gray contour lines of the arrival front,
and puts ESRI World Imagery underneath via contextily.

Needs numpy, matplotlib, contextily and xyzservices. On this machine the `geo`
conda environment has them:

    /Users/karl/anaconda3/envs/geo/bin/python scripts/plot_arrival.py \
        --case examples/Dish/0_Dish/output/dish -o dish_arrival.png

Design notes
------------
Arrival time is a *magnitude* field, so the fill is a sequential ramp, not a
categorical one. A heat ramp (`inferno`) is used deliberately: multi-hue
sequential is normally wrong, but "semantic heat" is one of the two sanctioned
exceptions, and it is paired with a colorbar as that exception requires. Rainbow
ramps (jet, turbo, rainbow, hsv) are refused -- they are not perceptually
uniform and invent structure that is not in the data.

Contours are gray on purpose. The fill already owns the colour channel for
magnitude; the contours add a second, redundant reading of the same field, so
they stay neutral instead of competing. Their gray ramps light-to-dark with
time, which is the "gradient in gray" that makes the front direction readable
without a second legend.
"""

import argparse
import os
import sys

try:
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import cm, colors
    from matplotlib.patheffects import withStroke
except ImportError as e:
    sys.exit("error: %s\nThis script needs numpy and matplotlib. Try:\n"
             "  /Users/karl/anaconda3/envs/geo/bin/python %s ..."
             % (e, sys.argv[0]))

NODATA_FALLBACK = -9999.0
RAINBOW = {"jet", "turbo", "rainbow", "hsv", "gist_rainbow", "nipy_spectral"}

# Neutral ink, so text never wears the data colour.
INK = "#f5f5f3"
INK_MUTED = "#c3c2b7"
HALO = "#0b0b0b"


def read_ascii_grid(path):
    """Parse an ESRI ASCII grid. Header keys are matched case-insensitively
    because FARSITE writes NCOLS/NROWS while other tools write ncols/nrows."""
    with open(path) as fh:
        tokens = fh.read().split()
    hdr, i = {}, 0
    while i + 1 < len(tokens):
        key = tokens[i].lower()
        if key not in ("ncols", "nrows", "xllcorner", "yllcorner", "cellsize",
                       "nodata_value", "xllcenter", "yllcenter"):
            break
        hdr[key] = tokens[i + 1]
        i += 2
    for req in ("ncols", "nrows", "cellsize"):
        if req not in hdr:
            sys.exit("error: %s is missing '%s'; is it an ESRI ASCII grid?" % (path, req))
    ncols, nrows = int(hdr["ncols"]), int(hdr["nrows"])
    cell = float(hdr["cellsize"])
    if "xllcorner" in hdr:
        xll, yll = float(hdr["xllcorner"]), float(hdr["yllcorner"])
    else:                                    # centre-referenced variant
        xll = float(hdr["xllcenter"]) - cell / 2.0
        yll = float(hdr["yllcenter"]) - cell / 2.0
    nodata = float(hdr.get("nodata_value", NODATA_FALLBACK))

    vals = np.array(tokens[i:i + ncols * nrows], dtype=float)
    if vals.size != ncols * nrows:
        sys.exit("error: %s declares %dx%d = %d cells but holds %d values"
                 % (path, ncols, nrows, ncols * nrows, vals.size))
    # Row 0 of the file is the north edge.
    return vals.reshape(nrows, ncols), xll, yll, cell, nodata


def add_scalebar(ax, length_m, label, color=INK):
    """Simple cartographic scale bar in data (metre) units."""
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    pad_x = (x1 - x0) * 0.06
    pad_y = (y1 - y0) * 0.06
    xs = x0 + pad_x
    ys = y0 + pad_y
    ax.plot([xs, xs + length_m], [ys, ys], color=color, lw=3, solid_capstyle="butt",
            path_effects=[withStroke(linewidth=5, foreground=HALO)], zorder=6)
    ax.text(xs + length_m / 2.0, ys + pad_y * 0.22, label, color=color,
            ha="center", va="bottom", fontsize=8,
            path_effects=[withStroke(linewidth=2.5, foreground=HALO)], zorder=6)


def nice_scalebar_length(span_m):
    for cand in (100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000):
        if cand >= span_m * 0.12:
            return cand, ("%d m" % cand if cand < 1000 else "%g km" % (cand / 1000.0))
    return 100000, "100 km"


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--case", help="output base path, e.g. .../output/dish "
                                    "(appends _ArrivalTime.asc)")
    src.add_argument("--grid", help="path to an arrival-time .asc directly")
    p.add_argument("-o", "--output", default="arrival_time.png")
    p.add_argument("--crs", default="EPSG:5070",
                   help="CRS of the grid, for the basemap fetch (default: EPSG:5070)")
    p.add_argument("--units", choices=["minutes", "hours"], default="hours",
                   help="FARSITE writes minutes; 'hours' converts (default: hours)")
    p.add_argument("--cmap", default="inferno", help="sequential colormap (default: inferno)")
    p.add_argument("--cmap-floor", type=float, default=0.18,
                   help="skip this fraction of the colormap's dark end so the "
                        "earliest cells stay visible over imagery (default: 0.18)")
    p.add_argument("--alpha", type=float, default=0.80,
                   help="fill opacity over the imagery (default: 0.80)")
    p.add_argument("--contour-levels", type=int, default=6,
                   help="number of contour lines, 0 to disable (default: 10)")
    p.add_argument("--contour-labels", action="store_true", help="label the contours inline")
    p.add_argument("--max-labelled-levels", type=int, default=4,
                   help="at most this many contour levels get inline labels "
                        "(default: 4)")
    p.add_argument("--no-basemap", action="store_true", help="skip the imagery fetch")
    p.add_argument("--zoom", default="auto", help="basemap zoom (default: auto)")
    p.add_argument("--attribution", default="",
                   help="basemap attribution text. Default empty, matching the "
                        "project's existing usage; note ESRI's terms of use if "
                        "you publish the figure.")
    p.add_argument("--title", help="plot title (default: derived from the case name)")
    p.add_argument("--perimeters", help="optional perimeter shapefile to outline")
    p.add_argument("--dpi", type=int, default=200)
    p.add_argument("--figsize", nargs=2, type=float, default=[9.0, 8.0])
    args = p.parse_args()

    if args.cmap.rstrip("_r") in RAINBOW:
        sys.exit("error: %r is a rainbow colormap. Magnitude needs a perceptually\n"
                 "uniform ramp -- try inferno, magma, plasma, viridis or cividis."
                 % args.cmap)

    grid_path = args.grid or ("%s_ArrivalTime.asc" % args.case)
    if not os.path.exists(grid_path):
        sys.exit("error: no arrival-time grid at %s" % grid_path)

    arr, xll, yll, cell, nodata = read_ascii_grid(grid_path)
    nrows, ncols = arr.shape

    # Unburned and nodata cells drop out so the imagery shows through.
    # Note the threshold is < 0, not <= 0: FARSITE writes arrival time 0 for the
    # cells of the initial ignition perimeter, so masking zero punches a hole
    # exactly where the fire started.
    burned = np.ma.masked_where((arr == nodata) | (arr < 0), arr)
    if burned.count() == 0:
        sys.exit("error: %s has no burned cells to plot" % grid_path)
    if args.units == "hours":
        burned = burned / 60.0
        unit_label = "hours since ignition"
    else:
        unit_label = "minutes since ignition"

    # Cell edges for pcolormesh; y descends from the north edge.
    x_edges = xll + np.arange(ncols + 1) * cell
    y_edges = (yll + nrows * cell) - np.arange(nrows + 1) * cell

    print("grid:    %s" % grid_path)
    print("  %d x %d cells at %g m, %d burned" % (ncols, nrows, cell, burned.count()))
    print("  extent  x[%.1f, %.1f] y[%.1f, %.1f] (%s)"
          % (x_edges[0], x_edges[-1], y_edges[-1], y_edges[0], args.crs))
    print("  arrival %.2f to %.2f %s" % (burned.min(), burned.max(), args.units))

    # Crop the view to the burned area plus a margin, rather than the whole
    # landscape -- a 15 km landscape with a 2 km fire is mostly empty map.
    rows, cols = np.where(~burned.mask)
    margin = max(10, int(0.25 * max(np.ptp(rows) + 1, np.ptp(cols) + 1)))
    r0, r1 = max(0, rows.min() - margin), min(nrows, rows.max() + margin + 1)
    c0, c1 = max(0, cols.min() - margin), min(ncols, cols.max() + margin + 1)
    view = (x_edges[c0], x_edges[c1], y_edges[r1], y_edges[r0])

    fig, ax = plt.subplots(figsize=tuple(args.figsize))
    fig.patch.set_facecolor("#111110")
    ax.set_facecolor("#111110")

    if not args.no_basemap:
        try:
            import contextily as ctx
            import xyzservices.providers as xyz
            ax.set_xlim(view[0], view[1])
            ax.set_ylim(view[2], view[3])
            # Create ESRI basemap (World Imagery)
            tile_provider = xyz.Esri.WorldImagery
            ctx.add_basemap(ax, source=tile_provider, zoom=args.zoom, crs=args.crs,
                            attribution=args.attribution)
            print("basemap: Esri.WorldImagery at zoom %s" % args.zoom)
        except ImportError as e:
            print("warning: no basemap (%s). Install contextily + xyzservices, or "
                  "pass --no-basemap." % e)
        except Exception as e:
            print("warning: basemap fetch failed (%s); continuing without it." % e)

    vmin, vmax = float(burned.min()), float(burned.max())
    norm = colors.Normalize(vmin=vmin, vmax=vmax)
    # Heat ramps bottom out at near-black, which disappears against dark
    # satellite imagery -- the earliest-arrival core reads as a hole in the map.
    # Truncating the low end keeps it a visible dark violet instead.
    base_cmap = (matplotlib.colormaps[args.cmap] if hasattr(matplotlib, "colormaps")
                 else cm.get_cmap(args.cmap))
    if args.cmap_floor > 0.0:
        fill_cmap = colors.LinearSegmentedColormap.from_list(
            "%s_clipped" % args.cmap,
            base_cmap(np.linspace(args.cmap_floor, 1.0, 256)))
    else:
        fill_cmap = base_cmap
    mesh = ax.pcolormesh(x_edges, y_edges, burned, cmap=fill_cmap, norm=norm,
                         alpha=args.alpha, shading="flat", zorder=3,
                         edgecolors="none", rasterized=True)

    # Contours of the same field: a redundant, neutral reading of the front.
    if args.contour_levels > 0:
        xc = 0.5 * (x_edges[:-1] + x_edges[1:])
        yc = 0.5 * (y_edges[:-1] + y_edges[1:])
        levels = np.linspace(vmin, vmax, args.contour_levels + 2)[1:-1]
        # Gray ramp, clipped away from both ends: pure white and pure black both
        # disappear against satellite imagery.
        gray_ramp = (matplotlib.colormaps["Greys_r"] if hasattr(matplotlib, "colormaps")
                     else cm.get_cmap("Greys_r"))
        grays = gray_ramp(np.linspace(0.30, 1.0, len(levels)))
        filled = burned.filled(np.nan)
        cs = ax.contour(xc, yc, filled, levels=levels, colors=grays,
                        linewidths=1.1, zorder=4)
        halo = [withStroke(linewidth=2.0, foreground=HALO, alpha=0.55)]
        if hasattr(cs, "collections"):          # matplotlib < 3.10
            for coll in cs.collections:
                coll.set_path_effects(halo)
        else:                                    # ContourSet is itself a Collection
            cs.set_path_effects(halo)
        if args.contour_labels:
            # Label a subset: a convoluted front relabels the same level many
            # times over and the numbers end up on top of each other.
            step = max(1, len(levels) // args.max_labelled_levels)
            lbl = ax.clabel(cs, levels=levels[::step], inline=True, fontsize=7,
                            inline_spacing=6,
                            fmt=("%.1f" if args.units == "hours" else "%.0f"))
            for t in lbl:
                t.set_color(INK)
                t.set_path_effects([withStroke(linewidth=2.0, foreground=HALO)])
        print("contours: %d levels, %.2f to %.2f" % (len(levels), levels[0], levels[-1]))

    if args.perimeters and os.path.exists(args.perimeters):
        try:
            from osgeo import ogr
            ogr.UseExceptions()
            ds = ogr.Open(args.perimeters)
            n = 0
            for feat in ds.GetLayer(0):
                g = feat.GetGeometryRef()
                for gi in range(g.GetGeometryCount() or 1):
                    ring = g.GetGeometryRef(gi) if g.GetGeometryCount() else g
                    if ring is None:
                        continue
                    pts = np.array([[ring.GetX(k), ring.GetY(k)]
                                    for k in range(ring.GetPointCount())])
                    if len(pts) > 1:
                        ax.plot(pts[:, 0], pts[:, 1], color=INK_MUTED, lw=0.7,
                                alpha=0.85, zorder=5)
                        n += 1
            ds = None
            print("perimeters: %d rings from %s" % (n, os.path.basename(args.perimeters)))
        except Exception as e:
            print("warning: could not draw perimeters (%s)" % e)

    ax.set_xlim(view[0], view[1])
    ax.set_ylim(view[2], view[3])
    ax.set_aspect("equal")

    # Recessive chrome. Axis numbers in kilometres: raw EPSG:5070 metres are
    # seven-digit values that crowd the axis and tell the reader nothing.
    ax.tick_params(colors=INK_MUTED, labelsize=8, length=3, width=0.6)
    for s in ax.spines.values():
        s.set_color("#3a3a38")
        s.set_linewidth(0.8)
    # Pick the decimal count from the span: rounding a 4 km extent to whole
    # kilometres yields repeated tick labels like "1912, 1912, 1912".
    span_km = max(view[1] - view[0], view[3] - view[2]) / 1000.0
    dec = 0 if span_km >= 40 else (1 if span_km >= 4 else 2)
    fmt = matplotlib.ticker.FuncFormatter(lambda v, _, d=dec: "%.*f" % (d, v / 1000.0))
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=6, steps=[1, 2, 5, 10]))
        axis.set_major_formatter(fmt)
    for t in ax.get_xticklabels():
        t.set_rotation(0)
    ax.set_xlabel("%s easting (km)" % args.crs, color=INK_MUTED, fontsize=8)
    ax.set_ylabel("%s northing (km)" % args.crs, color=INK_MUTED, fontsize=8)

    name = args.title or ("Fire arrival time — %s"
                          % os.path.basename(args.case or grid_path).replace("_ArrivalTime.asc", ""))
    ax.set_title(name, color=INK, fontsize=13, pad=12, loc="left")

    length, lbl = nice_scalebar_length(view[1] - view[0])
    add_scalebar(ax, length, lbl)

    cbar = fig.colorbar(mesh, ax=ax, fraction=0.040, pad=0.02, extend="neither")
    cbar.set_label(unit_label, color=INK_MUTED, fontsize=9)
    cbar.ax.yaxis.set_tick_params(color=INK_MUTED, labelsize=8)
    cbar.outline.set_edgecolor("#3a3a38")
    plt.setp(plt.getp(cbar.ax.axes, "yticklabels"), color=INK_MUTED)
    cbar.solids.set_alpha(1.0)

    fig.tight_layout()
    fig.savefig(args.output, dpi=args.dpi, facecolor=fig.get_facecolor())
    print("\nwrote %s" % args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
