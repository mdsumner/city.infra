# sfnetworks and the suite (wkpool, cdtr, traipse)

Status: investigation note, 2026-11-07. Nothing filed upstream.

Sources read: sfnetworks `main` at fe6edf8 (2026-05-14, version 0.9.9.9000,
the v1.0 rewrite), its 27 open issues (titles and opening posts; comment
threads were not reachable from the sandbox), and wkpool, cdtr and traipse
`main`. The Hobart/Tasmania numbers come from
[network_checks.py](network_checks.py) run on the tas-transport-web layer
(snapshot 2025-07-10). R was not available in the sandbox, so no sfnetworks
or wkpool code was actually run; the checks reproduce the same logic in
Python.

## 1. How sfnetworks works

An `sfnetwork` is a tidygraph `tbl_graph` (igraph underneath) whose node
table is an `sf` POINT table and whose edge table carries integer
`from`/`to` columns and, optionally, `sf` LINESTRING geometry ("explicit"
edges). `validate_network()` insists that node coordinates equal edge
boundary coordinates exactly, same CRS, same precision.

**Construction from lines** (`as_sfnetwork.sf` -> `create_from_spatial_lines`,
`R/create.R`):

1. `sfheaders::sf_to_df()` explodes every edge to one row per coordinate.
2. First and last coordinate of each linestring are the candidate nodes.
3. Node identity is decided by `st_match_points_df()` (`R/utils.R`): each
   coordinate column is `round()`ed to `precision_digits()` (12 decimal
   digits when the sf precision is 0), the columns are `paste()`d into a
   string key, and `match(key, unique(key))` gives the node index.
4. Optional `subdivide = TRUE` uses the same string keys over *all*
   coordinates and splits an edge wherever an interior vertex shares a key
   with any other vertex.

Only LINESTRING is accepted. MULTILINESTRING has to be cast first, which
changes the row count.

**Cleaning** is done with "spatial morphers" (`R/morphers.R`):

| Morpher | What it does | Implementation |
|---|---|---|
| `to_spatial_subdivision` | split edges at interior vertices shared with other edges | string keys again (`R/subdivide.R`); does not node at bare crossings (#134) |
| `to_spatial_smooth` | remove pseudo nodes (degree 2, or in=out=1 when directed), merge their edges | `igraph::decompose()` of the pseudo-node subgraph, then one R closure plus `sf::st_line_merge()` per set (`R/smooth.R`); `protect` and `require_equal` keep nodes |
| `to_spatial_contracted` | merge groups of nodes (e.g. dbscan clusters) into one | `R/contract.R` |
| `to_spatial_simple` | drop loops and multi-edges | `R/simplify.R` |
| `to_spatial_unique` | merge nodes at the same location | |
| `st_network_blend` | snap external points onto edges, inserting nodes | `R/blend.R` |
| `st_network_join` | join two networks on node location | `R/join.R` |

There are converters both ways for dodgr (`R/dodgr.R`), spatstat linnet/psp
and tidygraph.

## 2. Project health (worth knowing before offering anything)

- **#294 (2026-08-10, CRAN):** the maintainer address bounced and CRAN asked
  urgently for a new version. I could not reach CRAN from the sandbox to see
  whether the package has since been archived *(to verify)*.
- **#292:** igraph heads-up that `all_shortest_paths()` will validate a scalar
  `from`; `R/paths.R` needs a change.
- **#287 (high priority, v1.0):** `to_spatial_subdivision()` assigns wrong
  node indices in some cases (Cape Town workshop data). Still open.
- `main` has been in the v1.0 rewrite since late 2024; last commit May 2026.

So the project is active in spirit but thinly maintained. Small, low-risk,
well-tested help is likely welcome; large new dependencies probably are not.

## 3. The Tasmanian network through sfnetworks' eyes

All 284,831 LIST transport segments (Road 270,628, Track 12,110, Rail 2,083,
Route 10), 4.37 M vertices, 252,727 nodes from `node0`/`node1`. Whole run is
about 15 s in Python.

| Check | Result | sfnetworks equivalent |
|---|---|---|
| Node ids vs exact coordinate identity of endpoints | 252,727 = 252,727, one to one | `as_sfnetwork()` node matching |
| String keys at 12 / 7 / 6 decimal digits (lon/lat) | 252,727 / 252,689 / 252,680 nodes | `precision` argument |
| Distinct nodes closer than 1 cm / 10 cm / 1 m | 48 / 53 / 67 pairs | `to_spatial_unique`, contraction |
| Self-loops | 809 | `to_spatial_simple` |
| Multi-edges | 2,331 node pairs, 2,336 extra edges | `to_spatial_simple` |
| Degree 1 / 2 / 3 / 4 / 5+ | 77,342 / 47,116 / 116,158 / 11,357 / 754 | |
| Connected components | 244; largest has 95.8% of nodes; 139 are a single edge | `group_components()` |
| Edges after smoothing all pseudo nodes | about 237,700 (from 284,831) | `to_spatial_smooth()` |
| ... keeping nodes where `TRANS_TYPE` changes | about 238,900 | `require_equal = "TRANS_TYPE"` |
| ... keeping nodes where `PRI_NAME` changes | about 241,800 | `require_equal = "PRI_NAME"` |
| Interior vertices exactly on another segment's node | 259, on 222 segments | `to_spatial_subdivision()` |
| Dead ends within 2 m of another segment's vertex | 297 | undershoots (#115) |
| Segment pairs crossing with no shared vertex | 742 (385 road-road, 194 road-track, 112 road-rail) | noding at crossings (#134) |
| MULTILINESTRING rows with more than one part | 54 | would need casting before `as_sfnetwork()` |

Reading: the LIST layer is well noded already. The node ids are exactly
"same coordinate, same node", which is what `wkpool::merge_coincident(tolerance = 0)`
gives. The real cleaning questions are small and local: about 50 sub-centimetre
near misses, 259 vertex T-junctions, about 300 possible undershoots, and 742
crossings that are mostly grade separations and must *not* be noded blindly.
Those are a good review list for the mapper, not a job for a global tolerance.

## 4. Where the suite fits

Ranked by value to sfnetworks and by how easy it is to land.

### 4.1 Exact node identity without string keys (wkpool idea, #130)

sfnetworks decides "same node" by formatting rounded doubles into strings.
That is the root of #130 (coordinates that differ in the last bits get the
same key in one place and different keys in another, giving the igraph
"cannot add negative number of vertices" error), and it is slow on millions
of vertices. wkpool replaced exactly this pattern: `vertex_groups()` in
`R/merge_coincident.R` groups on the double bit pattern with
`vctrs::vec_group_id()`, and its comment notes that formatted-text keys
"silently collapsed doubles differing beyond ~15 significant digits".

Behaviour also differs by CRS today: 12 decimal digits plus the integer part
fits in 15 significant digits for lon/lat, but not for projected metres
(5,250,000.x), where `paste()` falls back to 15 significant digits.

Offer: a small PR to `st_match_points_df()` / `st_duplicated_points_df()`
using `vctrs::vec_group_id()` on the rounded numeric columns (or on raw
columns when precision is 0). vctrs is already a transitive dependency via
dplyr, so no new package is needed. This fixes a 2021 bug, is easy to test
and does not ask them to adopt wkpool.

### 4.2 Build from pre-noded ids, and a file format (#291)

`sfnetwork(nodes, edges, force = TRUE)` with integer `from`/`to` skips all
coordinate matching. Our web layer already carries `node0`/`node1`, so the
recipe is: read the GeoParquet, take node points from the first/last
coordinates keyed by id, rename `node0`/`node1` to `from`/`to`, and call
`sfnetwork()`. wkpool can produce the same columns for any line data
(`establish_topology()` + `merge_coincident()` gives `.vx0`/`.vx1` with
`.feature` provenance).

#291 asks for `write_sfnetwork()`/`read_sfnetwork()` "possibly as a single
.gpkg with two tables". tas-transport-web is a worked example of an
alternative: one GeoParquet edge table with integer endpoint ids, nodes
derivable, bbox covering for range reads, readable in the browser. That is
worth a comment on #291 with a link once the layer has a README.

### 4.3 Smoothing and subdivision as one arc walk (wkpool)

`wkpool::find_arcs()` walks maximal chains through degree-2 vertices in C++.
On a pool built from *all* segment vertices (not just endpoints), any
interior vertex shared with another segment gets the same `.vx` and so has
degree above 2. That means one `merge_coincident()` + `find_nodes()` +
`find_arcs()` pass yields the subdivided *and* smoothed network that
sfnetworks reaches with two morphers, each doing per-coordinate string keys,
and smoothing doing one `st_line_merge()` per pseudo-node set.

Gaps before this could back `to_spatial_smooth()`:

- `protect` and `require_equal`: arcs are attribute-blind. Those vertices
  need to be forced to nodes. `find_arcs()` currently takes only `x` and
  `quotient`, so this would be a new argument (for example `nodes =`).
- Directed networks: sfnetworks' directed pseudo node is in = out = 1; arcs
  are undirected.
- Provenance: sfnetworks can keep `.orig_data`; wkpool keeps `.feature` per
  segment, which is enough to rebuild it.

This is the strongest technical fit but the biggest ask upstream. A
benchmark on the Tasmanian network (sfnetworks morphers vs the wkpool walk,
same output count) is the evidence to bring first; #6 "Benchmarking" is the
natural place.

### 4.4 Snapping and undershoots (wkpool, with care, #115)

`merge_coincident(tolerance > 0)` snaps to a grid (`round(v / tol) * tol`).
That is fast and deterministic, but two points 1 mm apart can still fall in
different cells, and chains of near points do not merge transitively. For
the 48 sub-centimetre pairs in the LIST data a radius pass (kd-tree pairs
then union-find) is the right tool. wkpool's "observe, don't correct"
principle fits #115 well: report candidates (the 297 dead ends near another
segment) and let the user decide, rather than silently snapping.

### 4.5 Noding at bare crossings (cdtr, #134)

`cdtr::cdt_triangulate(..., intersect = "resolve", edge_sources = TRUE)`
inserts crossing points between constraint segments and reports, for every
output edge, which input segment it came from and in which direction. That
is a noding engine with provenance, which is exactly what #134 asks for
without going through GEOS `st_node()` (which loses attributes). Caveat
from the Tasmanian check: of 742 crossings most are probably bridges and
overpasses, so any "node all crossings" option must be opt-in and ideally
filterable by attribute (rail vs road, a bridge flag).

### 4.6 Turn angles (traipse, #290)

#290 asks about turn restrictions. traipse computes geodesic bearings and
turn angles along tracks (`track_bearing()`, `track_turn()`). Applied to
consecutive edge pairs at each node it gives a turn-angle table that both
sfnetworks and dodgr turn penalties could use. Lower priority; the dodgr/m4ra
thread is the better home for this.

## 5. Suggested next steps

1. Run the R side for real: build an `sfnetwork` from tas-transport-web with
   `from`/`to` (4.2), time `to_spatial_smooth()` against `wkpool::find_arcs()`
   on Hobart LGA and on all of Tasmania (4.3). Needs an R environment with
   CRAN access.
2. Draft the #130 PR (4.1) with a test using the coordinates from the issue.
   Ask Michael before opening anything upstream.
3. Turn the review list in section 3 (near misses, T-junctions, undershoots,
   crossings) into a mapper layer so they can be checked against the
   orthophoto.
4. Check sfnetworks' CRAN status (#294) before investing upstream; if it is
   archived, a working, faster sfnetworks matters more to the R spatial
   network community, and also to UrbanAnalyst users who meet it via dodgr's
   converters.
