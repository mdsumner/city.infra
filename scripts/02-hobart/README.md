# Hobart pipeline scripts (plan step 2)

Walk and walk + public transport time from every Greater Hobart address to
the Royal Hobart Hospital. Write-up and friction log:
[docs/02-hobart-pipeline.md](../../docs/02-hobart-pipeline.md).

## Inputs

| File | Source |
|---|---|
| `segs.parquet` | `https://data.source.coop/mdsumner/tas-transport-web/2025-07-10/list_transport_segments.parquet` |
| `addr.parquet` | `https://data.source.coop/mdsumner/thelist-opendata/geoparquet/list_address_points/statewide/list_address_points_statewide.parquet` |
| `gtfs.zip` | Tasmanian statewide GTFS (Public Transport Tasmania), from the archive in [`scripts/gtfs/`](../gtfs/) (`tas-statewide/<snapshot>/gtfs.zip`). The November 2026 run used `tas-statewide/2026-11-08` (version 25); the first run used version 20 from a third-party copy. |

## Run

```sh
python3 01_extract_network.py segs.parquet edges.csv     # needs pyarrow
python3 02_extract_addresses.py addr.parquet addresses.csv
Rscript 03_access_rhh.R edges.csv addresses.csv gtfs.zip out/
```

The analysis day is the first Tuesday on or after the later of today and the
feed's start date with no `calendar_dates` removals; set `GTFS_DATE=YYYYMMDD`
to fix it (the first run used `GTFS_DATE=20260804`).

For the published map with geographic context (land and water, LGA
boundaries, suburb names, overview inset, scale), fetch the theLIST context
layers once and run `04_map_access.R` on the outputs:

```sh
B=https://data.source.coop/mdsumner/thelist-opendata/geoparquet
mkdir -p context
curl -o context/lga.parquet $B/list_local_govt_areas/statewide/list_local_govt_areas_statewide.parquet
for l in hobart glenorchy clarence kingborough brighton derwent_valley huon_valley sorell; do
  curl -o context/loc_$l.parquet $B/list_locality_postcode/$l/list_locality_postcode_$l.parquet
done
Rscript 04_map_access.R out/rhh-access-addresses.csv.gz context/ out/run.json out/rhh-access-map.png
```

`04_map_access.R` needs arrow, sf, data.table and jsonlite. Travel time uses a
single-hue blue ramp (light is quick, dark is slow) checked as an ordinal
palette; water is grey so it is not read as data.

`03_access_rhh.R` needs dodgr, gtfsrouter, m4ra, data.table, sf and
geodist. It takes about 30 seconds on 4 cores. Set `USE_M4RA=1` to build the
stop-to-stop matrix with `m4ra_gtfs_traveltimes()` instead of
`gtfsrouter::gtfs_traveltimes()` (faster, but capped at 60 minutes).

Outputs in `out/`: `rhh-access-by-lga.csv`, `rhh-access-by-locality.csv`,
`rhh-access-by-stop.csv` (best time from each boarding stop to the hospital
door), `rhh-access-map.png`, `run.json` (feed hash, date, parameters), and `rhh-access-addresses.csv.gz` (one row per
address; not committed).
