# Compare wkpool::find_faces() with neighbourhoods::network_cycles() on the
# same graph: Greater Hobart LIST transport segments as straight edges
# between LIST node ids (the form network_cycles() works on).
#
# Both get identical input: largest connected component, no self-loops,
# one edge per node pair. Cycles are compared as sets of edge ids.
#
# Usage: Rscript compare_network_cycles.R edges.csv outdir
#   edges.csv from scripts/02-hobart/01_extract_network.py

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

# ---- neighbourhoods::network_cycles ------------------------------------------
x <- data.frame (
    .vx0 = as.character (e$node0), .vx1 = as.character (e$node1),
    .vx0_x = e$x0, .vx0_y = e$y0, .vx1_x = e$x1, .vx1_y = e$y1,
    edge_ = as.character (e$TRANSEG_ID),
    d = e$COMP_LEN
)
t_nc <- system.time (nc <- network_cycles (x))
# network_cycles() reports reversed edges as "<id>_rev"
key_nc <- vapply (nc, function (p) paste (sort (unique (sub ("_rev$", "", p$edge_))), collapse = ","), character (1))
message (sprintf ("network_cycles: %d cycles in %.1f s", length (nc), t_nc [["elapsed"]]))

# ---- wkpool::find_faces ------------------------------------------------------
v <- unique (rbind (e [, .(id = node0, x = x0, y = y0)], e [, .(id = node1, x = x1, y = y1)]), by = "id")
v [, .vx := seq_len (.N)]
pool <- new_wkpool (
    vertices = data.frame (.vx = v$.vx, x = v$x, y = v$y),
    vx0 = v$.vx [match (e$node0, v$id)], vx1 = v$.vx [match (e$node1, v$id)]
)
t_ff <- system.time (ff <- find_faces (pool))
ff <- as.data.table (ff)
ff [, edge_ := as.character (e$TRANSEG_ID [.segment])]
bounded <- ff [.face_area > 0 & !.bridge]
key_ff <- bounded [, .(key = paste (sort (unique (edge_)), collapse = ",")), by = .face]
message (sprintf ("find_faces: %d bounded faces in %.2f s", nrow (key_ff), t_ff [["elapsed"]]))

# ---- compare -------------------------------------------------------------------
key_outer <- ff [.face_area < 0 & !.bridge, .(key = paste (sort (unique (edge_)), collapse = ",")), by = .face]$key
face_sets <- strsplit (key_ff$key, ",")
edge_faces <- split (rep (seq_along (face_sets), lengths (face_sets)), unlist (face_sets))
# is an extra network_cycles() cycle exactly the outline of a group of faces?
is_union <- function (k) {
    ce <- strsplit (k, ",") [[1]]
    cand <- unique (unlist (edge_faces [intersect (ce, names (edge_faces))]))
    if (length (cand) < 2) return (FALSE)
    inside <- cand [vapply (face_sets [cand], function (fs) mean (fs %in% ce) > 0.5, logical (1))]
    if (length (inside) < 2) return (FALSE)
    tab <- table (unlist (face_sets [inside]))
    setequal (names (tab) [tab == 1], ce)
}
# is a find_faces() face self-crossing (a sign it lies where chords cross)?
vxy <- data.table (.vx = v$.vx, x = v$x, y = v$y)
self_crossing <- function (fid) {
    # dead ends are walked out and back, so leave them out of the outline
    w <- ff [.face == fid & !.bridge]
    xy <- vxy [match (c (w$.vx0, w$.vx0 [1]), .vx), cbind (x, y)]
    !sf::st_is_simple (sf::st_sfc (sf::st_linestring (xy)))
}

both <- intersect (key_nc, key_ff$key)
only_nc <- setdiff (key_nc, key_ff$key)
only_ff <- setdiff (key_ff$key, key_nc)
dup_nc <- sum (duplicated (key_nc))
nv <- nrow (v)
ne <- nrow (e)
nf <- length (unique (ff$.face))

only_nc_union <- vapply (only_nc, is_union, logical (1))
ff_only_ids <- key_ff$.face [key_ff$key %in% only_ff]
ff_only_crossing <- vapply (ff_only_ids, self_crossing, logical (1))
all_crossing <- vapply (key_ff$.face, self_crossing, logical (1))
summary <- list (
    edges = ne, vertices = nv, parallel_dropped = n_parallel,
    euler_V_minus_E_plus_F = nv - ne + nf, genus = (2 - (nv - ne + nf)) / 2,
    network_cycles = length (key_nc), network_cycles_duplicates = dup_nc,
    network_cycles_seconds = t_nc [["elapsed"]],
    find_faces_bounded = nrow (key_ff), find_faces_seconds = t_ff [["elapsed"]],
    in_both = length (both), only_network_cycles = length (only_nc),
    only_network_cycles_outer = sum (only_nc %in% key_outer),
    only_network_cycles_union_of_faces = sum (only_nc_union),
    only_find_faces = length (only_ff),
    only_find_faces_self_crossing = sum (ff_only_crossing),
    find_faces_self_crossing_all = sum (all_crossing)
)
jsonlite::write_json (summary, file.path (outdir, "compare-summary.json"), auto_unbox = TRUE, pretty = TRUE)
fwrite (data.table (method = "network_cycles", edges = only_nc), file.path (outdir, "only-network-cycles.csv"))
fwrite (data.table (method = "find_faces", edges = only_ff), file.path (outdir, "only-find-faces.csv"))
str (summary)
