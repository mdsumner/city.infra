# Step 1: the UrbanAnalyst stack, and where our tools plug in

Status: first draft, 2026-11-06. Based on reading the default branch of each
UrbanAnalyst repo (DESCRIPTION, README, and selected source files). Points
marked *(to verify)* are inferences that need checking with the maintainers
or by running code.

## 1. The pipeline in one paragraph

UA builds a city in five stages. Global inputs come first: OpenStreetMap
extracted with `osmium-tool`, GHSL population density rasters, and NASA
elevation. Per-city inputs are a GTFS feed and a set of polygons carrying a
socio-demographic variable (unemployment or a richer disadvantage index);
those polygons also define the city's extent. OSM is converted to
`osmdata_sc` (silicate `SC` form) and weighted into a routable graph by
**dodgr**. **gtfsrouter** handles the transit leg, and **m4ra** combines
walk / bike / transit into many-to-many travel times and ratios against
private car. **uaengine** computes the per-vertex variables (transport
times, bicycle infrastructure share, natural space share, parking, school
distance), weights by population density, and aggregates into the social
polygons. The website then uses **uamutations** (Rust) to "mutate" one city's
variable distribution toward another's, and **docs** (Sphinx) explains it all.

```
OSM (osmium) -> osmdata_sc -> dodgr graph --+
GTFS ------------------------> gtfsrouter --+--> m4ra --> uaengine --> polygons --> website
GHSL raster, elevation, social polygons ----------------^                         uamutations
```

## 2. Repo by repo

| Repo | Role | Lang | Main inputs | Status badge | Open issues |
|------|------|------|-------------|--------------|-------------|
| dodgr | Weighted directed graphs, many-to-many distances, flows | R + C++ (Rcpp, RcppParallel) | `sf`, `sc`/`SC` street nets; `highway`-style type column + weighting profile | active | 18 |
| gtfsrouter | GTFS routing, timetables, isochrones | R + C++ (Rcpp) | GTFS zip | active | 15 |
| m4ra | Multi-modal many-to-many travel times, ratios to car | R + C++ | `osmdata_sc` + GTFS | concept | 2 |
| uaengine | Variables + aggregation for the UA site | R + C++ (cpp11) | m4ra output, GHSL GeoTIFF (via `raster`), social polygons (`sf`); needs `osmium-tool` | concept | 2 |
| ttcalib | Calibrate car times (signals, turns) to empirical data | R | Uber Movement travel times (Santiago) | concept | 1 |
| neighbourhoods | Minimal cycles of a planar street graph; LTN identification and scoring | R + C++ (cpp11) | contracted, undirected dodgr graph; popdens file | concept | 4 |
| superblock | Superblock potential: space, parking, open space | R + C++ (Rcpp) | OSM via `sb_osmdata_extract()` | (none shown) | 4 |
| naturoute | Routing through natural spaces | R + C++ (cpp11) | `sf` | active | 0 |
| spatialcluster | Spatially contiguous clustering (REDCAP) | R + C++ (RcppArmadillo) | points + covariance/distance; optional boundary via **cdtr** | wip | 8 |
| us-hardship-index | US hardship index from ACS | R | US Census API | active | 0 |
| uamutations | Mutate city A toward city B's variable distribution | Rust | per-city JSON from uaengine | (CI only) | 1 |
| docs | Platform documentation | Python (Sphinx, rst) | - | - | 3 |

Issue counts are from the GitHub search API on 2026-11-06.

## 3. Language boundaries

- **R is the integration layer.** Every analytic component is an R package;
  the website is fed by R-produced files.
- **C++ does the heavy lifting**, split between two binding styles: Rcpp
  (+ RcppParallel / RcppThread) in the older core (dodgr, gtfsrouter, m4ra,
  superblock, spatialcluster) and cpp11 in the newer pieces (uaengine,
  neighbourhoods, naturoute). wkpool and cdtr are both cpp11, so they sit
  naturally with the newer half.
- **Rust** appears only in uamutations, which runs on the site side and reads
  JSON written by uaengine.
- **System tools**: `osmium-tool` is a hard requirement of uaengine.
- **Geometry layer**: OSM arrives in silicate `SC` form (`osmdata_sc`), which
  is close kin to wkpool's vertex-pool model. Rasters go through the retired
  `raster` package; polygons through `sf` / `sfheaders`.

## 4. Where our tools plug in

### wkpool

- **theLIST network -> dodgr graph.** The LIST transport segments are already
  cleanly noded linestrings. `establish_topology()` + `merge_coincident()`
  gives a vertex pool and segment table, which is exactly dodgr's `from` /
  `to` / `d` edge form, with `.feature` carrying the LIST segment id. This
  would be the first non-OSM route into the UA stack.
- **Cycles for neighbourhoods.** `neighbourhoods::network_cycles()` finds
  minimal cycles with its own C++ walk. Checked 2026-11-08: wkpool's
  `find_cycles()` recovers polygon rings, it does not find the faces of a
  line network, so the contribution is to add face enumeration to wkpool
  alongside its existing `find_neighbours()`, `find_shared_edges()` and
  `classify_cycles()`. Draft: [wkpool-neighbourhoods.md](wkpool-neighbourhoods.md).
- **Adjacency for aggregation.** Edge and vertex neighbour relations between
  social polygons or between neighbourhood cycles are what LTN aggregation and
  spatialcluster contiguity both need.

### cdtr

- **Already used.** spatialcluster suggests cdtr for a boundary-constrained
  neighbour graph (`R/edges.R`). Getting that onto CRAN-ready footing helps
  directly.
- **Population onto the network.** uaengine assigns GHSL cells to vertices via
  `raster::rasterToPoints()`. A triangulated mesh of the social polygons with
  `cdt_interpolate()` is an option for areal-to-network weighting *(to
  verify)*, and the `raster` dependency itself is a candidate to move to
  `terra` or GDAL-direct reads.

### traipse

- **Geodesic metrics.** dodgr uses geodist for edge lengths; traipse adds
  bearings and turn angles, relevant to dodgr's turn penalties.
- **Empirical calibration.** ttcalib calibrates against Uber Movement, which
  is no longer available *(to verify)*. Local GPS traces (rides, walks,
  bus AVL) processed with traipse are a possible replacement source, and
  Hobart could be the test case.

### theLIST open data (hypertidy/thelist-opendata)

- Authoritative transport network, cadastre, land tenure and parcels for
  Tasmania. These could substitute for, or validate, OSM in Hobart.
- The LIST road classes would need a dodgr `wt_profile_file` that maps them
  to the `highway`-type categories dodgr weights on.

### city.infra cycle mapper

- Produces human-curated groupings of LIST segments with descriptions: a
  ground-truth layer for "what is this infrastructure actually like".
- uaengine's bicycle-infrastructure variable relies on OSM tags. Curated
  Hobart groupings could validate or replace it locally, and the same
  pattern extends to accessibility (kerb ramps, footpath condition, crossing
  quality) where OSM is thin.
- Groupings of segments are also natural candidate inputs for `superblock`
  and `neighbourhoods` (pick a patch, score it).

## 5. Gaps noticed while reading

Framing and docs:

- Several core repos carry a "concept" status badge (m4ra, uaengine, ttcalib,
  neighbourhoods), so a newcomer can't tell what is production.
- `neighbourhoods` still points badges at `ATFutures-labs` and has
  `Remotes: ATFutures/dodgr`.
- Small README nits (for example `install.packges` in spatialcluster, a stale
  version string in gtfsrouter).
- There is no single "how the pieces fit" page; section 1 above is a first cut.

Wiring:

- No documented path for non-OSM networks into m4ra / uaengine.
- Social variable is US-specific (`us-hardship-index`). The Australian
  analogue is ABS SEIFA (IRSD) at SA1 or SA2, which would make any
  Australian city possible.
- `raster` dependency in uaengine.

Data:

- Calibration data source (Uber Movement) appears gone *(to verify)*.
- GTFS for Hobart: Metro Tasmania feed availability and licence *(to
  verify)*.

## 6. Questions to take to the UA maintainers

1. Which repos are on the critical path for the live site, and which are
   exploratory?
2. Is a non-OSM network input (theLIST) welcome in dodgr / m4ra, or better
   kept as a conversion layer outside?
3. Would wkpool as a shared topology dependency for `neighbourhoods` be
   welcome, or is the in-house C++ preferred?
4. What replaced Uber Movement for calibration, if anything?
5. Is there appetite for an Australian city (Hobart) on the platform, and what
   is the minimum data bar?
