#!/usr/bin/env python3
"""Compare two ESRI shapefiles (.shp) numerically.

Byte comparison is too strict: the coordinates FARSITE writes are doubles, so a
last-bit floating-point difference changes the file while the geometry is
effectively identical. This compares structure exactly (record count, record
numbers, shape types, part/point counts) and coordinates with a relative
tolerance.

Exit 0 if equivalent, 1 otherwise. Usage: compare_shp.py A.shp B.shp [rtol]
"""
import struct
import sys

# Shape types that store a bbox + parts before their point array.
MULTI = {3: "PolyLine", 5: "Polygon", 23: "PolyLineM", 25: "PolygonM"}
POINT = {1: "Point", 11: "PointZ", 21: "PointM"}


def records(buf):
    """Yield (recno, shapetype, [parts], [(x, y), ...]) for each record."""
    pos = 100
    n = len(buf)
    while pos + 8 <= n:
        recno, wlen = struct.unpack_from(">ii", buf, pos)
        body, end = pos + 8, pos + 8 + 2 * wlen
        if end > n:
            raise ValueError("truncated record %d" % recno)
        (stype,) = struct.unpack_from("<i", buf, body)
        pts, parts = [], []
        if stype in POINT:
            pts = [struct.unpack_from("<2d", buf, body + 4)]
        elif stype in MULTI:
            nparts, npoints = struct.unpack_from("<ii", buf, body + 4 + 32)
            off = body + 4 + 32 + 8
            parts = list(struct.unpack_from("<%di" % nparts, buf, off))
            off += 4 * nparts
            pts = [struct.unpack_from("<2d", buf, off + 16 * i) for i in range(npoints)]
        elif stype != 0:  # 0 = null shape
            raise ValueError("unsupported shape type %d" % stype)
        yield recno, stype, parts, pts
        pos = end


def close(a, b, rtol):
    return abs(a - b) <= rtol * max(abs(a), abs(b), 1.0)


def main():
    pa, pb = sys.argv[1], sys.argv[2]
    rtol = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-6
    A, B = open(pa, "rb").read(), open(pb, "rb").read()

    for name, off, fmt in (("file code", 0, ">i"), ("version", 28, "<i"), ("shape type", 32, "<i")):
        va = struct.unpack_from(fmt, A, off)[0]
        vb = struct.unpack_from(fmt, B, off)[0]
        if va != vb:
            print("%s differs: %s vs %s" % (name, va, vb))
            return 1

    worst = 0.0
    for i, (ba, bb) in enumerate(zip(struct.unpack_from("<4d", A, 36), struct.unpack_from("<4d", B, 36))):
        if not close(ba, bb, rtol):
            print("header bbox[%d] differs: %.17g vs %.17g" % (i, ba, bb))
            return 1
        worst = max(worst, abs(ba - bb) / max(abs(bb), 1.0))

    try:
        ra, rb = list(records(A)), list(records(B))
    except ValueError as e:
        print("parse error: %s" % e)
        return 1

    if len(ra) != len(rb):
        print("record count differs: %d vs %d" % (len(ra), len(rb)))
        return 1

    ndiff = 0
    for (na, ta, pa_, xa), (nb, tb, pb_, xb) in zip(ra, rb):
        if (na, ta, pa_) != (nb, tb, pb_):
            print("record %d structure differs (no/type/parts)" % na)
            return 1
        if len(xa) != len(xb):
            print("record %d point count differs: %d vs %d" % (na, len(xa), len(xb)))
            return 1
        for (x1, y1), (x2, y2) in zip(xa, xb):
            for u, v in ((x1, x2), (y1, y2)):
                if u != v:
                    ndiff += 1
                    rel = abs(u - v) / max(abs(v), 1.0)
                    worst = max(worst, rel)
                    if not close(u, v, rtol):
                        print("record %d coord differs beyond tol: %.17g vs %.17g" % (na, u, v))
                        return 1

    if ndiff:
        print("%d coords differ, max rel %.3g (within tol %.3g)" % (ndiff, worst, rtol))
    return 0


if __name__ == "__main__":
    sys.exit(main())
