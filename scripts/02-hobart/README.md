# Hobart pipeline scripts (plan step 2)

Walk and walk + public transport time from every Greater Hobart address to
the Royal Hobart Hospital. Write-up and friction log:
[docs/02-hobart-pipeline.md](../../docs/02-hobart-pipeline.md).

## Inputs

| File | Source |
|---|---|
| `segs.parquet` | `https://data.source.coop/mdsumner/tas-transport-web/2025-07-10/list_transport_segments.parquet` |
| `addr.parquet` | `https://data.source.coop/mdsumner/thelist-opendata/geoparquet/list_address_points/statewide/list_address_points_statewide.parquet` |
| `gtfs.zip` | Tasmanian statewide GTFS (Public Transport Tasmania). The run here used feed version 20, valid 2026-06-03 to 2026-08-31 (see the doc for where it came from). |

## Run

```sh
python3 01_extract_network.py segs.parquet edges.csv     # needs pyarrow
python3 02_extract_addresses.py addr.parquet addresses.csv
Rscript 03_access_rhh.R edges.csv addresses.csv gtfs.zip out/
```

`03_access_rhh.R` needs dodgr, gtfsrouter, m4ra, data.table, sf and
geodist. It takes about 30 seconds on 4 cores. Set `USE_M4RA=1` to build the
stop-to-stop matrix with `m4ra_gtfs_traveltimes()` instead of
`gtfsrouter::gtfs_traveltimes()` (faster, but capped at 60 minutes).

Outputs in `out/`: `rhh-access-by-lga.csv`, `rhh-access-by-locality.csv`,
`rhh-access-map.png`, and `rhh-access-addresses.csv.gz` (one row per
address; not committed).
