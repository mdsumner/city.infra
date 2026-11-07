# thelist-opendata: where it got to

Checked 2026-11-07 against the live store
(https://source.coop/mdsumner/thelist-opendata, listed anonymously via
`s3://us-west-2.opendata.source.coop/mdsumner/thelist-opendata/`) and the
GitHub repo hypertidy/thelist-opendata.

## Done

- Pipeline steps 00 to 03 have all run: crawl 2026-09-10, convert 2026-09-17,
  upload 2026-10-13.
- Store holds 3116 objects, 20 GiB: `geoparquet/` (219 products, 12.7 GiB),
  `cog/` (50 products, 7.3 GiB), `parquet/` (83), `manifest/`, `README.md`.
- `manifest/converted.csv` on the store: 3024 conversions, all `ok` (30 via
  the FlatGeobuf hop).
- The store README is thorough: layout, attribution wording, licence (CC BY
  3.0 AU), read examples for DuckDB, GDAL, R and Python, and the conversion
  caveats.
- **CORS works**: `data.source.coop` answers range requests with
  `access-control-allow-origin: *`, so a browser app can read the GeoParquet
  directly.
- `list_transport_segments/` is there for all 29 LGAs (69 MiB total) with a
  `list_transport_segments_union.vrt`. This can replace the splh.transport
  cache.

## Loose ends

1. **`raw/` is missing.** The store README says the source zips are alongside
   (183 GB), but there are no `raw/` objects. Either upload them, or change the
   README to say raw is not mirrored and point at the LIST URLs.
2. **Boundary files are duplicated.** The README says `municipality_<lga>` is
   converted once under `geoparquet/municipality/` (32 files there), but 1321
   copies also sit inside other product folders, each with its own
   `municipality_union.vrt`. Either remove those or fix the README.
3. **Stray object**: `README.md/README.md` (a second copy under a "folder"
   named README.md).
4. **14 sources not converted.** These are the expected ones: the CHM and
   Slope 2 m rasters, which need a local copy, and the known empty or odd
   zips, such as the CFEV Clarence and Derwent catchments, storm tide
   Kingborough, EPA documents, NCH blueberries and barley v2. Also
   `INTERIM_PLANNING_SCHEME_OVERLAY_STATEWIDE`, `SES_DESIGN_FLOOD_EVENT_STATEWIDE`,
   `LIST_SHELLMAP_HARVEST_AREA_BOUNDARIES_STATEWIDE` and
   `LIST_TREE_CANOPY_HEIGHT_STATEWIDE`, which are worth a look.
5. **GitHub repo lags the store.** `manifest/converted.csv` is not committed,
   and the repo README does not link the source.coop product or carry the
   reading examples.
6. **No refresh.** Nothing re-runs the crawl. The data is the 2026-09-10
   snapshot. A scheduled GitHub Action, using a source.coop service account,
   would need to run 00 to 03 with etag skipping.
7. **No snapshot versioning.** source.coop has none. If ids such as
   `TRANSEG_ID` matter for long-lived references (the mapper groupings),
   decide between dated prefixes and recording `converted.csv` etags with
   each grouping.
8. Still EPSG:28355 throughout, by design. The web layer for the mapper will
   be a derived product (EPSG:4326 plus node ids), kept separate from the
   faithful mirror.
