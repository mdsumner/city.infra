# Cycle mapper: review and a path to a no-R, Pages-hosted version

Status: draft, 2026-11-06. Based on a full read of `apps/cycle_mapper/app.R`
(about 1000 lines) and of `hypertidy/thelist-opendata`. The app was not run
here; points marked *(to verify)* need a check against the real data.

## 1. Short answer

Yes. Nothing in the app needs a server. It is a filter, a small graph walk
and a click-to-select map. Saving is a CSV. All of that runs in the
browser. The one heavy part is building endpoint adjacency with nested R
loops on every search, and that should be **precomputed once, offline**,
not ported. With topology precomputed, the browser side becomes small.

Recommended shape:

- **Data (offline, R + GDAL):** the network is published as a normal
  thelist-opendata product on S3, plus one derived web-ready artefact.
- **Site (Quarto, GitHub Pages):** a Quarto website for city.infra, with the
  mapper as one page built in Observable JS (OJS). Quarto supports OJS
  natively, and its built-in `DuckDBClient` can query Parquet over HTTP from S3.
- **Saving:** browser `localStorage` for in-progress work, plus CSV
  export and import. Crowdsourced groupings come back as CSVs, by PR or by
  upload, and merge on `TRANSEG_ID` as they do now.

## 2. Folding the cache into thelist-opendata

The cache at `mdsumner/splh.transport` releases
(`LIST_TRANSPORT_SEGMENTS.parquet`, 84 MB) duplicates a product that
thelist-opendata already maps. `manifest/layers.csv` lists
`LIST_TRANSPORT_SEGMENTS` (and `LIST_TRANSPORT_NODES`) per LGA, in EPSG:28355
(Hobart LGA has 9471 segments, 26 fields). `R/03_convert.R` already writes
per-LGA GeoParquet plus a per-product union to `DEST` (for example
`/vsis3/bucket/thelist`). So:

1. Let `03_convert.R` produce `geoparquet/list_transport_segments/...` and
   the union, like any other product. Retire the release asset.
2. Add a **derived** layer, kept separate from the faithful mirror, since
   the README's rule is "nothing is materialised that can be derived". The web app needs:
   - geometry in EPSG:4326 (so the browser never needs PROJ),
   - `TRANSEG_ID`, `PRI_NAME`, `SEC_NAME`, `COMP_LEN`,
   - `node0`, `node1`: integer endpoint node ids from
     `wkpool::merge_coincident()` (or from `LIST_TRANSPORT_NODES` if the
     segments carry node references *(to verify)*),
   - rows sorted spatially (Hilbert), with a GeoParquet 1.1 bbox covering
     column, so range reads over HTTP only fetch the area in view.
3. Optionally, a PMTiles vector-tile version of the same layer, so the whole
   network can be drawn at any zoom with MapLibre and clicked by id.
4. The bucket needs CORS enabled for the Pages origin.

**CRS (resolved):** the cache recipe (`splh.transport/data-raw/transport-segments.R`)
does no reprojection, so the cache is native EPSG:28355 (GDA94 / MGA55). The
app assumed EPSG:7855 (GDA2020), which shifted display by about 1.5 m; it now
uses 28355. See [data/README.md](../data/README.md).

## 3. Review of the current R app

What works well and should carry over:

- The street syntax (`Argyle [Federal, Burnett], Weld`) is a strong idea. It
  is human-readable, reproducible, and could itself be stored as the
  provenance of a grouping.
- Output is just `id, description, transeg_id`, which is deterministic and
  easy to merge.
- Context segments either side of a bounded section are a nice touch for
  people checking a selection on the map.

Issues, roughly by impact:

1. **Adjacency is O(n^2) in R loops and is recomputed per search**
   (`order_segments_by_network()` and `filter_between_cross_streets()`, which
   also builds it twice). It is slow on long streets, and the 1 m / 5 m
   tolerances are fragile. Precomputed node ids fix both problems, and the
   comparisons become integer equality.
2. **Range between cross streets relies on a greedy walk order.** When a
   street branches or has a gap, `order_segments_by_network()` jumps to the
   nearest unvisited segment, so "segments between index i and j" can
   include the wrong pieces. A shortest path on the street's own subgraph
   between the two intersection nodes is more robust and no harder in JS.
3. **No provenance on saved groupings.** There's no author, no date, no
   network version, and no stored query. `TRANSEG_ID` stability across LIST
   releases is the key data question for crowdsourcing *(to verify)*. If ids
   can change, each grouping should record the network snapshot it was made on.
4. **Shared CSV in the working directory.** On a multi-user Shiny server, users
   overwrite each other, and saving an existing id silently replaces it.
   (This goes away in the static version, where each user saves locally.)
5. **Search with an empty query draws every segment in the 10 km box**, one
   `addPolylines()` call each. That will be slow. Vector tiles or a single
   GeoJSON layer avoid it.
6. **Table and map selection fight.** The table observer replaces
   `selected()` with the table's rows, so segments picked on the map are
   lost after a click in the table *(to verify by running)*.
7. Small points:
   - The "10 km" Hobart box uses 1/111 km per degree for longitude too. At
     -42.9 latitude it is about 13.7 km east-west.
   - The comment documents `[A:]` / `[:B]`, but the parser and help text
     use `[A,]`. A single `[A]` falls through to "whole street" without a
     message.
   - The app is titled "cycle infrastructure", but nothing in it is cycle
     specific. A neutral name (segment grouping for any infrastructure
     quality) widens it to accessibility, footpaths and crossings.

## 4. Proposed next steps

1. Add `LIST_TRANSPORT_SEGMENTS` conversion and the derived web layer
   (4326, node ids, Hilbert-sorted) in thelist-opendata, to S3 with CORS.
2. Set up a Quarto website in city.infra (`_quarto.yml`, landing page, the
   existing docs as pages) published to GitHub Pages by a GitHub Action.
3. Port the mapper to a Quarto OJS page: DuckDB-WASM for name and bbox
   queries, Leaflet or MapLibre for the map, then the syntax parser and
   graph walk in JS, with localStorage and CSV import and export.
4. Keep the R app working against the new S3 layer in the meantime, since
   that is a one-line change of `NETWORK_URL` plus the CRS fix.
