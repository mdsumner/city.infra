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
9. **Water crossings are missing from the `tas-transport-web` segment
   geometry** (found 2026-11-08, web layer snapshot `2025-07-10`). Across
   the state, 609 segments have at least 50 m less geometry than their
   official `COMP_LEN`, and 54 are MULTILINESTRINGs whose parts have gaps.
   In Greater Hobart all three road crossings of the River Derwent are
   affected:

   | TRANSEG_ID | Road | `COMP_LEN` (m) | geometry (m) | crossing |
   |---|---|---:|---:|---|
   | 4298506 | Tasman Highway | 1,264 | 266 | Tasman Bridge |
   | 4298327, 4298328 | Goodwood Road | 993 | 131, 133 | Bowen Bridge |
   | 6410527, 6411293 | Midland Highway | 1,283, 1,292 | 971, 976 | Bridgewater |

   Statewide the list includes ferries, the Batman Highway (Batman Bridge)
   and the Tasman Highway at Sorell. Likely cause *(to verify against the
   source zips)*: LIST is distributed per LGA, the parts of a segment over
   water lie in no LGA, and they are lost when the per-LGA files are
   combined; the `lga` column then lists only the land LGAs (the Tasman
   Bridge segment is tagged `clarence` only). Effects: on the full geometry
   the Greater Hobart network splits into eastern and western shores; graphs
   built from end-node ids join the gaps with straight edges but the Tasman
   Bridge edge does not reach the Hobart side
   ([docs/02-hobart-pipeline.md](../docs/02-hobart-pipeline.md), friction
   point 9; [docs/wkpool-neighbourhoods.md](../docs/wkpool-neighbourhoods.md)).
   Fix in the web-layer build: take each `TRANSEG_ID`'s full geometry once
   from the statewide source rather than combining per-LGA clips, and add a
   check of geometry length against `COMP_LEN`.
