# Data used by city.infra

## Web layer used by the segment mapper

The browser mapper (`mapper.qmd`, `site/mapper/`) reads
[tas-transport-web](https://source.coop/mdsumner/tas-transport-web):
`latest.json` names the current dated snapshot, which holds
`list_transport_segments.parquet` (EPSG:4326, one row per `TRANSEG_ID`,
`node0`/`node1` endpoint node ids, `snapshot` column, GeoParquet 1.1 bbox
covering), a PMTiles copy, and `snapshot.json` recording the source zip etags.
Saved groupings record the snapshot they were made on.

## Transport segments network (Shiny app)

| | |
|---|---|
| File | `list_transport_segments_statewide.parquet` (55 MB, GeoParquet 1.1) |
| URL | https://data.source.coop/mdsumner/thelist-opendata/geoparquet/list_transport_segments/statewide/list_transport_segments_statewide.parquet |
| S3 | `s3://us-west-2.opendata.source.coop/mdsumner/thelist-opendata/geoparquet/list_transport_segments/statewide/` (anonymous, path-style) |
| Built | 2026-11-20, from the per-LGA conversion of the 2026-09-10 crawl |
| Provenance | `manifest/converted.csv` on the store (row `lga-union:list_transport_segments/...`); per-LGA source etags in the same file |
| Source | LIST open data, `LIST_TRANSPORT_SEGMENTS_<LGA>.zip`, all 29 LGAs |
| Method | [hypertidy/thelist-opendata](https://github.com/hypertidy/thelist-opendata): per-LGA FileGDB to GeoParquet, then pieces clipped at LGA boundaries rejoined to one row per `TRANSEG_ID` (284831 rows from 286157 pieces); column `lga` lists the source LGAs |
| CRS | Native, not reprojected: EPSG:28355 (GDA94 / MGA zone 55) |
| Key fields used | `TRANSEG_ID`, `PRI_NAME`, `SEC_NAME`, `COMP_LEN`, geometry `geometry` (WKB) |
| Licence | "LIST Transport Segments" from theLIST (c) State of Tasmania, CC BY 3.0 AU |

This replaces the interim cache
`https://github.com/mdsumner/splh.transport/releases/download/latest/LIST_TRANSPORT_SEGMENTS.parquet`
(release 2025-08-29), which had the same 284831 `TRANSEG_ID`s but kept the
LGA-clipped pieces as separate rows (286157 rows). Groupings that reference
`TRANSEG_ID`s carry over unchanged.

Known gap, inherited from the per-LGA source: parts of a segment that lie
over water fall in no LGA file, so they are not in the statewide file either.
14 segments are more than 50 m shorter than their `COMP_LEN`, including the
three Derwent road crossings (Tasman Bridge `4298506`, Bowen Bridge
`4298327`/`4298328`, Bridgewater `6410527`/`6411293`), the Batman Bridge
(`6111927`), five `TRANS_TYPE = Route` segments (likely ferries) and
three short track bridges. See
[notes/thelist-opendata-status.md](../notes/thelist-opendata-status.md).

The store has no snapshot versioning: the file is replaced in place when the
pipeline re-runs. Record the build date (or the etag from
`manifest/converted.csv`) alongside any shared groupings.
