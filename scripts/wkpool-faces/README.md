# wkpool find_faces() against neighbourhoods::network_cycles()

Write-up: [docs/wkpool-neighbourhoods.md](../../docs/wkpool-neighbourhoods.md).

| Script | What |
|---|---|
| `compare_planar_grids.R` | Jittered grids (planar by construction) where the true face count is known; writes `planar-grids.csv`. |
| `compare_network_cycles.R` | Greater Hobart LIST network as edges between node ids; writes `compare-summary.json` and the cycles found by only one method. |

Needs wkpool with `find_faces()` (hypertidy/wkpool#14), neighbourhoods
(install from GitHub; it needs caret), igraph, data.table and sf.

```sh
Rscript compare_planar_grids.R out/
Rscript compare_network_cycles.R ../02-hobart/edges.csv out/   # edges.csv from 01_extract_network.py
```
