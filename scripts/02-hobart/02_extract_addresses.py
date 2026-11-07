"""Extract LIST address points as a flat CSV (EPSG:28355 easting/northing).

Usage: python3 02_extract_addresses.py ADDRESS_POINTS.parquet OUT.csv
ADDRESS_POINTS.parquet is
https://data.source.coop/mdsumner/thelist-opendata/geoparquet/list_address_points/statewide/list_address_points_statewide.parquet
"""
import sys

import pyarrow.csv as pcsv
import pyarrow.parquet as pq

COLS = ["GEOCODE_ID", "GEO_TYPE", "EASTING", "NORTHING", "SITE_NAME",
        "ST_NO_FROM", "STREET", "ST_TY_CODE", "LOCALITY", "POSTCODE", "RURAL_ADD"]

t = pq.read_table(sys.argv[1], columns=COLS)
pcsv.write_csv(t, sys.argv[2])
print(f"{t.num_rows} address points -> {sys.argv[2]}")
