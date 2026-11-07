# Data used by city.infra

## Web layer used by the segment mapper

The browser mapper (`mapper.qmd`, `site/mapper/`) reads
[tas-transport-web](https://source.coop/mdsumner/tas-transport-web):
`latest.json` names the current dated snapshot, which holds
`list_transport_segments.parquet` (EPSG:4326, one row per `TRANSEG_ID`,
`node0`/`node1` endpoint node ids, `snapshot` column, GeoParquet 1.1 bbox
covering), a PMTiles copy, and `snapshot.json` recording the source zip etags.
Saved groupings record the snapshot they were made on.

## Transport segments network cache (Shiny app, interim)

| | |
|---|---|
| File | `LIST_TRANSPORT_SEGMENTS.parquet` (84 MB, ZSTD) |
| Current URL | https://github.com/mdsumner/splh.transport/releases/download/latest/LIST_TRANSPORT_SEGMENTS.parquet |
| Last modified | 2025-08-29 (HTTP Last-Modified of the release asset) |
| Recipe | [splh.transport/data-raw/transport-segments.R](https://github.com/mdsumner/splh.transport/blob/main/data-raw/transport-segments.R) |
| Source | LIST open data, `LIST_TRANSPORT_SEGMENTS_<LGA>.zip`, all 29 LGAs, FileGDB layer `list_transport_segments_<lga>` |
| Method | `gdal vector convert` per LGA (streamed via `/vsizip//vsicurl/`), then `gdal vector concat --mode single` |
| CRS | Native, not reprojected: EPSG:28355 (GDA94 / MGA zone 55) per the thelist-opendata manifest |
| Key fields used | `TRANSEG_ID`, `PRI_NAME`, `SEC_NAME`, `COMP_LEN`, geometry `SHAPE` (WKB) |
| Licence | LIST open data, CC BY (attribution to Land Tasmania / theLIST) |

Planned home: a `LIST_TRANSPORT_SEGMENTS` product in
[hypertidy/thelist-opendata](https://github.com/hypertidy/thelist-opendata),
published on source.coop, plus a derived web layer for the static mapper (see
[docs/cycle-mapper-review.md](../docs/cycle-mapper-review.md)). When it moves,
update this file and `NETWORK_URL` in `apps/cycle_mapper/app.R`.

Groupings saved by the mapper (`cycle_mappings.csv`) reference `TRANSEG_ID`s
from this snapshot; record the snapshot date alongside any shared groupings.
