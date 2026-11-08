# Compare wkpool::find_faces() with neighbourhoods::network_cycles() on the
# full LIST geometry: every vertex of every Greater Hobart transport
# segment, so edges are the short straight pieces of each road line and
# the graph follows the real road shapes (no straight chords between
# junctions).
#
# Both get identical input: vertices shared by exact coordinates
# (wkpool::merge_coincident), largest connected component, no self-loops,
# one edge per vertex pair. Cycles are compared as sets of edge ids.
#
# Usage: Rscript compare_full_geometry.R segs.parquet outdir
#   segs.parquet: list_transport_segments.parquet from tas-transport-web
#   (see scripts/02-hobart/README.md)

source (file.path (dirname (sub ("^--file=", "", grep ("^--file=", commandArgs (FALSE), value = TRUE))), "compare_helpers.R"))
suppressPackageStartupMessages ({
    library (data.table)
    library (wkpool)
    library (neighbourhoods)
})

args <- commandArgs (trailingOnly = TRUE)
f_segs <- args [1]
outdir <- args [2]
dir.create (outdir, showWarnings = FALSE, recursive = TRUE)

lgas <- "hobart|glenorchy|clarence|kingborough|brighton"
not_walkable <- c ("Railway", "Railway Siding", "Ferry", "Tramway")
t <- arrow::read_parquet (f_segs, col_select = c ("TRANS_TYPE", "TRAN_CLASS", "lga", "geometry"))
keep <- grepl (lgas, t$lga) & t$TRANS_TYPE != "Rail" & !t$TRAN_CLASS %in% not_walkable
pool <- merge_coincident (establish_topology (wk::wkb (t$geometry [keep])))
s <- as.data.table (pool_segments (pool)) [.vx0 != .vx1]
s [, `:=` (lo = pmin (.vx0, .vx1), hi = pmax (.vx0, .vx1))]
n_parallel <- sum (duplicated (s [, .(lo, hi)]))
s <- s [!duplicated (s [, .(lo, hi)])]

g <- igraph::graph_from_data_frame (s [, .(as.character (.vx0), as.character (.vx1))], directed = FALSE)
cm <- igraph::components (g)
big <- as.integer (names (cm$membership) [cm$membership == which.max (cm$csize)])
s <- s [.vx0 %in% big]

pv <- pool_vertices (pool)
i0 <- match (s$.vx0, pv$.vx)
i1 <- match (s$.vx1, pv$.vx)
e <- data.table (
    node0 = s$.vx0, node1 = s$.vx1,
    x0 = pv$x [i0], y0 = pv$y [i0], x1 = pv$x [i1], y1 = pv$y [i1],
    edge_id = sprintf ("s%07d", seq_len (nrow (s)))
)
message (sprintf ("%d edges in the largest of %d components (%d parallel duplicates dropped)",
    nrow (e), cm$no, n_parallel))
compare_cycles (e, outdir, prefix = "full-geometry-",
    extra = list (input = "full LIST geometry", components_before = cm$no, parallel_dropped = n_parallel))
