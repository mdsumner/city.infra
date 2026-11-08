---
title: "Draft: wkpool contributing to neighbourhoods"
date: 2026-11-08
---

Status: **draft for discussion, nothing filed upstream.** Update
2026-11-08: `find_faces()` is built (hypertidy/wkpool#14) and compared with
`network_cycles()`; see "Results" below. Read against
`UrbanAnalyst/neighbourhoods` at `8e08513` (2026-11-06) and `hypertidy/wkpool`
at `ccbbeb9` (2026-10-25). Items marked *(to verify)* need a test before they
go anywhere.

## Summary

- **A correction first.** The stack map ([01-ua-stack.md](01-ua-stack.md))
  says wkpool "already has `find_cycles()`" to compare with
  `network_cycles()`. That is wrong as written: wkpool's `find_cycles()`
  recovers the rings of polygon input from path provenance. It does not find
  the faces (enclosed blocks) of a line network, which is what `neighbourhoods`
  needs. The contribution is to **add** face finding to wkpool, not to swap in
  something that exists.
- **The main offer:** a deterministic planar face enumeration in wkpool
  (half-edge walk with a rotation system), which would replace the
  multi-pass, hash-and-deduplicate search in `network_cycles()` with a single
  pass in which every directed edge belongs to exactly one face.
- **Two quick wins** that need no new dependency: a likely one-line bug in
  restart-edge naming, and replacing the quadratic `adjacent_cycles()` loop
  with a shared-edge index (wkpool's `find_neighbours()` pattern).

## What neighbourhoods does now

`neighbourhoods()` runs a nine-step pipeline on an OSM network in silicate
`SC` form:

1. dodgr: weight for motor cars, keep the largest component, contract
   degree-2 vertices, compute centrality, merge to an undirected graph.
2. `network_cycles()`: the minimal cycles (blocks) of that graph.
3. `adjacent_cycles()`: pairs of cycles that share edges.
4. `nbs_add_data()`: area, population density (from a GeoTIFF via
   raster / stars), centrality in and out, highway classes.
5. Filters (no primary / trunk / service shared edges, shared edge lower than
   its neighbours) and an estimated effect score; `ltn_score()` /
   `cut_nbs()` then test cuts by recomputing centrality.

How `network_cycles()` finds cycles (R/cycles.R, src/cycles.cpp,
src/clockwise.cpp):

- Duplicates every edge as `<id>_rev`, removes terminal vertices, then
  starts a left-turn trace from every edge, and a right-turn trace, and
  deduplicates the results by an MD5 hash of the sorted edge ids.
- Removes "isolated polygons" (groups joined to the rest by a single
  connection), because they act as "attractors" for the left-trace and stop
  it closing larger enclosing cycles. Then runs a second pass from "restart
  edges" (edges in only one cycle so far) and a final `cpp_reduce_paths()`.
- Vertex and edge ids are strings in `std::unordered_map`s throughout.

This works, but the extra passes, hashing and attractor handling are
symptoms of starting walks from edges and hoping to cover every face. The
standard construction covers every face by construction.

## What wkpool has, and does not have

Has, all on integer vertex ids into a shared pool:

- `establish_topology()` / `merge_coincident()`: exact vertex sharing
  (coincident coordinates become one vertex, no rounding to strings; the
  same approach as our fix for sfnetworks #130).
- Directed segments (`.vx0`, `.vx1`) with feature provenance.
- `vertex_degree()`, `find_nodes()`, `find_arcs()`: arc-node decomposition,
  which is degree-2 contraction (what `dodgr_contract_graph()` does).
- `shared_edge_index()`, `find_shared_edges()`, `find_neighbours()`:
  adjacency between features from shared edges or vertices, in one pass.
- `cycle_signed_area()`, `classify_cycles()`: winding and outer / hole
  roles.
- `tri_dissolve()`: a topological ring walk over a triangle mesh (turning
  through the triangle fan, no angle comparisons).

Does not have: face enumeration of an arbitrary segment graph. That is the
piece to build.

## Contributions, smallest first

### 1. Restart edges never include reverse edges *(to verify)*

`get_restart_edges()` (R/cycles.R:185) builds
`paste0 ("rev_", edges)`, but every reverse edge in the package is named
with a suffix, `paste0 (edge_, "_rev")` (R/network-preprocess.R:25, and
the C++ `cut_terminal_rev()`). So the reverse half of the restart set never
matches `x$edge_` and is dropped, and the second pass only restarts from
forward edges. Fix: `paste0 (edges, "_rev")`.

Before filing: a small network where a face is only reachable from a
reverse restart edge, run before and after.

### 1b. Short edge ids crash `network_cycles()`

`utils::cut_terminal_rev()` (src/utils.cpp) takes `size () - 4` on every
edge id, which underflows for ids shorter than four characters
(`basic_string::substr: __pos (which is 18446744073709551613)`). OSM ids
are long, so this never shows in normal use, but any other network with
short ids ("1", "2", ...) fails. Found while building the planar test
graphs below; fix: check the length before comparing the suffix.

### 2. `adjacent_cycles()` in one pass

Currently, for each cycle it scans every other cycle's edges with `%in%`,
so the cost grows with the square of the number of cycles. The same result
comes from one grouping: stack all (edge, cycle) pairs, keep edges that
appear in two cycles, and emit the pairs. This is exactly
`shared_edge_index()` / `find_neighbours()` in wkpool, and it can be written
inline (base R or data.table) without taking a wkpool dependency.

Also: `adjacent_cycles()` is exported and calls `dplyr::`, but dplyr is
only in Suggests, so it fails when dplyr is not installed. The rewrite drops
dplyr from this path.

### 3. Face enumeration in wkpool (the main piece)

A new wkpool function (working name `pool_faces()`):

1. Treat each undirected segment as two directed half-edges.
2. At every vertex, sort its outgoing half-edges by angle (the rotation
   system). One sort per vertex, by degree.
3. For a half-edge arriving at vertex *v*, the next half-edge on the same
   face is the one after its twin in *v*'s rotation order (always turn the
   same way).
4. Walk until every half-edge is used. Each half-edge lies on exactly one
   face, so there is no start-edge choice, no hashing and no deduplication.
5. Use the signed area (already in wkpool) to mark the outer face, and give
   each face its half-edges, vertices and the input segment ids.

What it would replace in `network_cycles()`: both trace directions, the
hash deduplication, isolated-polygon removal, restart edges and path
reduction. "Attractor" polygons are ordinary faces here; a part hanging off
by one connection, or a dead-end street, simply appears twice in the
boundary of the face that surrounds it and is trimmed from that face's ring
afterwards.

Pieces wkpool already supplies around it: exact vertex sharing before the
walk, `find_arcs()` for contraction, `shared_edge_index()` for face
adjacency (item 2 becomes a call), and `cycle_signed_area()` for winding.
The core walk would be C++ (cpp11, the same toolchain as
`neighbourhoods`).

Open design point: **street networks are not planar.** Bridges and
overpasses cross without a shared vertex. A rotation-system walk still
gives well-defined faces, but two faces can overlap where one road passes
over another. `network_cycles()` has the same property today. Options:
accept it (crossings without junctions are real), or offer a planarising
step that splits segments at crossings, with a flag on such split points.

### 4. Further out

- **Non-OSM networks.** wkpool can carry a LIST network (Hobart) straight
  into the graph form `neighbourhoods` needs, which is the same "bring your
  own network" gap as in m4ra and dodgr ([02-hobart-pipeline.md](02-hobart-pipeline.md),
  friction points 1 and 3).
- **Population onto faces.** `nbs_add_data()` reads a population-density
  GeoTIFF through raster and stars. With faces as polygons, census areal
  data (for Australia, the mesh block counts we already use) can be
  apportioned to faces by overlay instead, and the retired `raster`
  dependency could go. That is a separate conversation.

## Results

Harness: [scripts/wkpool-faces/](https://github.com/mdsumner/city.infra/tree/main/scripts/wkpool-faces);
result files in [wkpool-faces/](wkpool-faces/). Both methods get identical
input and cycles are compared as sets of edge ids. `find_faces()` faces are
the bounded ones (positive area), with dead-end edges left out.

**Planar test graphs.** Jittered grids (8, 20 and 40 vertices a side, 0, 15
and 30% of edges removed, 3 seeds each; up to 3,120 edges), where the true
number of faces is known from Euler's formula:

- `find_faces()` gives exactly the right faces in all 27 runs.
- `network_cycles()` finds every true face in 25 of 27 runs. It nearly
  always adds one extra cycle, which in 17 runs is exactly the outer
  boundary of the graph. In the other two runs it misses faces: in one, 2
  true faces are replaced by 2 cycles that each outline a group of faces; in
  the other, one face is missing.
- Time: `network_cycles()` up to 1.1 s on these graphs (`find_faces()` was
  not timed here; see Hobart below).

**Comparing like with like.** A face can have islands inside it: a loop at
the end of a dead-end road, a roundabout. Its walk goes round the islands'
outsides too, while `network_cycles()` returns the outer ring only. So each
`find_faces()` face is split into rings and only its outer ring (positive
area) is compared. Without this, faces with islands showed up as
"different" (97 faces on the full geometry).

**Greater Hobart, LIST, two versions of the same network.** Both methods get
identical input (largest connected component, no self-loops, one edge per
vertex pair).

- *Node graph*: straight edges between LIST node ids, the form
  `network_cycles()` is normally given (a contracted graph).
- *Full geometry*: every vertex of every road line, so edges follow the
  real road shapes.

| | node graph | full geometry |
|---|---:|---:|
| edges | 50,661 | 261,842 |
| genus (crossings without a shared vertex, roughly) | 133 | 37 |
| `network_cycles()` cycles | 6,436 (9 duplicates) | 4,145 (7 duplicates) |
| `find_faces()` bounded faces | 6,108 | 4,051 |
| found by both | 6,013 | 4,022 |
| only `network_cycles()` | 414 | 116 |
| of which outer boundary / outline of several faces | 14 / 22 | 1 / 14 |
| only `find_faces()` | 95 | 29 |
| of which self-crossing | 83 | 29 |
| `network_cycles()` time | about 20 s | 31 min |
| `find_faces()` time | 0.06 s | 0.4 s |

On the full geometry, 4,022 of 4,051 faces
(99.3%) are found by both, and **every face that
only `find_faces()` reports is self-crossing**: it lies where roads cross
without a junction (overpasses), the case that needs a decision either
way. 101
of the cycles only `network_cycles()` reports are not classified yet.
(The full-geometry `network_cycles()` result was rebuilt from the first
run's saved output to apply the outer-ring comparison, rather than
rerunning 31 minutes; the duplicate count is from the original run.)

**Caveat: the full geometry covers only the western shore.** Its largest
component has 261,842 of about 483,000 edges, because the three road
crossings of the River Derwent are broken in the `tas-transport-web` layer
(see [theLIST mirror status](../notes/thelist-opendata-status.md)). The node
graph is connected across the river because straight edges between node
ids jump the gaps. The comparison is still like for like (both methods get
the same component).

What this says for the proposal:

- On planar input `find_faces()` is exact where the current search is
  usually right but not always.
- It is about 300 times faster on the node graph (0.06 s against about
  20 s) and about 4,500 times faster on the full geometry (0.4 s against 31
  minutes), where the current search does not scale.
- Faces should be found on the full geometry and contracted afterwards (the
  face of each half-edge carries through contraction): the genus falls from
  133 to 37, so most crossing artefacts disappear, and the cost is trivial
  with `find_faces()`.
- The remaining real crossings need a decision (see the questions below).

## Next steps

- Classify the cycles only `network_cycles()` reports (101 on the full
  geometry), likely near the self-crossing faces.
- Rerun on the full geometry once the Derwent crossings are repaired in
  `tas-transport-web`.
- Draft issues for 1, 1b and 2, each with a minimal example (1b already has
  one in `compare_planar_grids.R`).
- A small adapter from `find_faces()` output to the `network_cycles()`
  result shape (a list of edge tables per cycle), so `adjacent_cycles()` and
  the rest of the pipeline can run on it unchanged for a side-by-side LTN
  result.

## Questions for the neighbourhoods maintainer

1. Is a wkpool dependency acceptable, or would you rather take the
   algorithm as C++ inside `neighbourhoods`?
2. Which network should the comparison use as the reference case?
3. How should overpasses be treated: as separate, overlapping faces, or
   split?
4. Is `neighbourhoods` still in active use for the platform, given its
   "concept" status badge and `ATFutures-labs` badge links?
