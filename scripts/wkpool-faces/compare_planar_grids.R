# Compare find_faces() and network_cycles() on graphs that are planar by
# construction (jittered grids with random edges removed), where the true
# number of bounded faces is known from Euler's formula: E - V + 1 for a
# connected planar graph.
#
# Usage: Rscript compare_planar_grids.R outdir

suppressPackageStartupMessages ({
    library (data.table)
    library (wkpool)
    library (neighbourhoods)
})
outdir <- commandArgs (trailingOnly = TRUE) [1]
dir.create (outdir, showWarnings = FALSE, recursive = TRUE)

grid_graph <- function (n, drop, seed) {
    set.seed (seed)
    g <- CJ (i = 0:(n - 1), j = 0:(n - 1))
    g [, `:=` (id = i + n * j + 1L, x = i + runif (.N, -0.3, 0.3), y = j + runif (.N, -0.3, 0.3))]
    h <- g [i < n - 1, .(a = id, b = id + 1L)]
    v <- g [j < n - 1, .(a = id, b = id + n)]
    e <- rbind (h, v)
    e <- e [runif (.N) > drop]
    # largest connected component only
    cm <- igraph::components (igraph::graph_from_data_frame (e [, .(as.character (a), as.character (b))], directed = FALSE))
    big <- as.integer (names (cm$membership) [cm$membership == which.max (cm$csize)])
    e <- e [a %in% big]
    # network_cycles() fails on edge ids shorter than 4 characters
    # (cut_terminal_rev() underflows), so pad them
    e [, edge := sprintf ("e%06d", seq_len (.N))]
    list (v = g, e = e)
}

one <- function (n, drop, seed) {
    gg <- grid_graph (n, drop, seed)
    v <- gg$v
    e <- gg$e
    x <- data.frame (
        .vx0 = as.character (e$a), .vx1 = as.character (e$b),
        .vx0_x = v$x [e$a], .vx0_y = v$y [e$a], .vx1_x = v$x [e$b], .vx1_y = v$y [e$b],
        edge_ = e$edge, d = 1
    )
    t_nc <- system.time (nc <- tryCatch (network_cycles (x), error = function (err) err))
    nc_ok <- !inherits (nc, "error")
    key_nc <- if (nc_ok) vapply (nc, function (p) paste (sort (unique (sub ("_rev$", "", p$edge_))), collapse = ","), character (1)) else character ()

    pool <- new_wkpool (vertices = data.frame (.vx = v$id, x = v$x, y = v$y), vx0 = e$a, vx1 = e$b)
    ff <- as.data.table (find_faces (pool))
    ff [, edge_ := e$edge [.segment]]
    key_ff <- ff [.face_area > 0 & !.bridge, .(key = paste (sort (unique (edge_)), collapse = ",")), by = .face]$key
    key_outer <- ff [.face_area < 0 & !.bridge, .(key = paste (sort (unique (edge_)), collapse = ",")), by = .face]$key
    only_nc_keys <- setdiff (key_nc, key_ff)
    # an extra cycle that is a union of true faces (its edges are the
    # symmetric difference of their boundaries) is a non-minimal cycle
    face_sets <- lapply (strsplit (key_ff, ","), identity)
    is_union <- vapply (only_nc_keys, function (k) {
        ce <- strsplit (k, ",") [[1]]
        inside <- vapply (face_sets, function (fs) all (fs %in% ce) || mean (fs %in% ce) > 0.5, logical (1))
        if (sum (inside) < 2) return (FALSE)
        tab <- table (unlist (face_sets [inside]))
        setequal (names (tab) [tab == 1], ce)
    }, logical (1))

    nv <- length (unique (c (e$a, e$b)))
    expected <- nrow (e) - nv + 1L
    data.table (
        n = n, drop = drop, seed = seed, edges = nrow (e), expected_faces = expected,
        find_faces = length (key_ff),
        network_cycles = if (nc_ok) length (key_nc) else NA_integer_,
        nc_duplicates = sum (duplicated (key_nc)),
        in_both = length (intersect (key_nc, key_ff)),
        only_nc = length (only_nc_keys),
        only_nc_outer = sum (only_nc_keys %in% key_outer),
        only_nc_union_of_faces = sum (is_union & !only_nc_keys %in% key_outer),
        only_ff = length (setdiff (key_ff, key_nc)),
        nc_seconds = t_nc [["elapsed"]],
        nc_error = if (nc_ok) "" else conditionMessage (nc)
    )
}

runs <- CJ (n = c (8L, 20L, 40L), drop = c (0, 0.15, 0.3), seed = 1:3)
res <- rbindlist (lapply (seq_len (nrow (runs)), function (k) one (runs$n [k], runs$drop [k], runs$seed [k])))
fwrite (res, file.path (outdir, "planar-grids.csv"))
print (res [, .(n, drop, seed, expected_faces, find_faces, network_cycles, in_both, only_nc, only_nc_outer, only_nc_union_of_faces, only_ff, nc_seconds)])
cat ("\nfind_faces == expected in", sum (res$find_faces == res$expected_faces), "of", nrow (res), "runs\n")
cat ("network_cycles finds every true face in", sum (res$only_ff == 0), "of", nrow (res), "runs\n")
cat ("network_cycles includes the outer boundary in", sum (res$only_nc_outer > 0), "of", nrow (res), "runs\n")
