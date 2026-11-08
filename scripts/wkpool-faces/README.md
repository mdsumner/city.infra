# wkpool find_faces() against neighbourhoods::network_cycles()

Write-up: [docs/wkpool-neighbourhoods.md](../../docs/wkpool-neighbourhoods.md).

| Script | What |
|---|---|
| `compare_planar_grids.R` | Jittered grids (planar by construction) where the true face count is known; writes `planar-grids.csv`. |
| `compare_network_cycles.R` | Greater Hobart LIST network as edges between node ids; writes `compare-summary.json` and the cycles found by only one method. |
| `compare_full_geometry.R` | The same network with every vertex of every road line (wkpool topology); writes `full-geometry-*` files. `network_cycles()` takes about 31 minutes here; its cycles are cached in `full-geometry-network-cycles-keys.csv` and reused if present. |
| `compare_helpers.R` | The shared comparison (outer rings of `find_faces()` faces against `network_cycles()` cycles, as edge-id sets). |

Needs wkpool with `find_faces()` (hypertidy/wkpool#14), neighbourhoods
(install from GitHub; it needs caret), igraph, data.table, sf and arrow.

```sh
Rscript compare_planar_grids.R out/
Rscript compare_network_cycles.R ../02-hobart/edges.csv out/   # edges.csv from 01_extract_network.py
Rscript compare_full_geometry.R segs.parquet out/              # list_transport_segments.parquet
```
