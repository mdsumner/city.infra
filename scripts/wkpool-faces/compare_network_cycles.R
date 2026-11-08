# Compare wkpool::find_faces() with neighbourhoods::network_cycles() on the
# same graph: Greater Hobart LIST transport segments as straight edges
# between LIST node ids (the form network_cycles() works on).
#
# Both get identical input: largest connected component, no self-loops,
# one edge per node pair. Cycles are compared as sets of edge ids.
#
# Usage: Rscript compare_network_cycles.R edges.csv outdir
#   edges.csv from scripts/02-hobart/01_extract_network.py

source (file.path (dirname (sub ("^--file=", "", grep ("^--file=", commandArgs (FALSE), value = TRUE))), "compare_helpers.R"))
suppressPackageStartupMessages ({
    library (data.table)
    library (wkpool)
    library (neighbourhoods)
})

args <- commandArgs (trailingOnly = TRUE)
f_edges <- args [1]
outdir <- args [2]
dir.create (outdir, showWarnings = FALSE, recursive = TRUE)

e <- fread (f_edges)
not_walkable <- c ("Railway", "Railway Siding", "Ferry", "Tramway")
e <- e [TRANS_TYPE != "Rail" & !TRAN_CLASS %in% not_walkable & node0 != node1]
e [, `:=` (lo = pmin (node0, node1), hi = pmax (node0, node1))]
n_parallel <- sum (duplicated (e [, .(lo, hi)]))
e <- e [!duplicated (e [, .(lo, hi)])]

# largest connected component
g <- igraph::graph_from_data_frame (e [, .(as.character (node0), as.character (node1))], directed = FALSE)
cm <- igraph::components (g)
big <- names (cm$membership) [cm$membership == which.max (cm$csize)]
e <- e [as.character (node0) %in% big]
message (sprintf ("%d edges in the largest component (%d parallel duplicates dropped)", nrow (e), n_parallel))
e [, edge_id := as.character (TRANSEG_ID)]
compare_cycles (e, outdir, extra = list (parallel_dropped = n_parallel))
