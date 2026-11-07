"""Extract Greater Hobart LIST transport segments into a flat edge table.

R on this machine has no Parquet reader (no arrow/nanoparquet, and the
system GDAL lacks the Parquet driver), so this step uses pyarrow and
writes a plain CSV that dodgr can take directly.

Usage: python3 01_extract_network.py SEGMENTS.parquet OUT.csv
SEGMENTS.parquet is list_transport_segments.parquet from
https://data.source.coop/mdsumner/tas-transport-web/<snapshot>/
"""
import struct
import sys

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.csv as pcsv
import pyarrow.parquet as pq

LGAS = {"hobart", "glenorchy", "clarence", "kingborough", "brighton"}


def ends(wkb):
    """First and last vertex of a WKB LineString or MultiLineString."""
    b = memoryview(wkb)
    le = "<" if b[0] == 1 else ">"
    gtype = struct.unpack(le + "I", b[1:5])[0] % 1000
    if gtype == 2:
        parts = [b]
    elif gtype == 5:
        n = struct.unpack(le + "I", b[5:9])[0]
        parts, off = [], 9
        for _ in range(n):
            ple = "<" if b[off] == 1 else ">"
            npt = struct.unpack(ple + "I", b[off + 5:off + 9])[0]
            parts.append(b[off:off + 9 + 16 * npt])
            off += 9 + 16 * npt
    else:
        raise ValueError(f"unexpected geometry type {gtype}")
    first, last = parts[0], parts[-1]
    fle = "<" if first[0] == 1 else ">"
    lle = "<" if last[0] == 1 else ">"
    nlast = struct.unpack(lle + "I", last[5:9])[0]
    x0, y0 = struct.unpack(fle + "dd", first[9:25])
    x1, y1 = struct.unpack(lle + "dd", last[9 + 16 * (nlast - 1):9 + 16 * nlast])
    return x0, y0, x1, y1


def main(src, out):
    t = pq.read_table(src, columns=["TRANSEG_ID", "PRI_NAME", "COMP_LEN", "TRANS_TYPE",
                                    "TRAN_CLASS", "lga", "node0", "node1", "snapshot", "geometry"])
    keep = [any(l in LGAS for l in (v or "").split(";")) for v in t["lga"].to_pylist()]
    t = t.filter(pa.array(keep))
    xy = [ends(g) for g in t["geometry"].to_pylist()]
    t = t.drop_columns(["geometry"])
    for i, name in enumerate(["x0", "y0", "x1", "y1"]):
        t = t.append_column(name, pa.array([p[i] for p in xy], pa.float64()))
    pcsv.write_csv(t, out)
    print(f"{t.num_rows} segments from {len(LGAS)} LGAs -> {out}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
