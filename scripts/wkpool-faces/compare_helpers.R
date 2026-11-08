# Shared comparison of neighbourhoods::network_cycles() and
# wkpool::find_faces() on one edge table. Sourced by the compare_*.R
# scripts.
#
# e: data.table with node0, node1 (vertex ids), x0, y0, x1, y1 and edge_id
#    (character, at least 4 characters: network_cycles() needs that), with
#    no self-loops, one edge per vertex pair, one connected component.
# Writes <prefix>compare-summary.json and the cycles found by one method.

compare_cycles <- function (e, outdir, prefix = "", extra = list ()) {

    # network_cycles() can take a long time (31 min on the full Hobart
    # geometry), so its cycles are cached by edge set and reused when the
    # cache exists
    f_cache <- file.path (outdir, paste0 (prefix, "network-cycles-keys.csv"))

    # ---- neighbourhoods::network_cycles ------------------------------------------
    x <- data.frame (
        .vx0 = as.character (e$node0), .vx1 = as.character (e$node1),
        .vx0_x = e$x0, .vx0_y = e$y0, .vx1_x = e$x1, .vx1_y = e$y1,
        edge_ = e$edge_id,
        d = 1
    )
    if (file.exists (f_cache)) {
        cache <- fread (f_cache)
        key_nc <- cache$edges
        t_nc <- c (elapsed = cache$seconds [1])
        message (sprintf ("network_cycles: %d cycles from cache (%.1f s when run)", length (key_nc), t_nc [["elapsed"]]))
    } else {
        t_nc <- system.time (nc <- network_cycles (x))
        # network_cycles() reports reversed edges as "<id>_rev"
        key_nc <- vapply (nc, function (p) paste (sort (unique (sub ("_rev$", "", p$edge_))), collapse = ","), character (1))
        fwrite (data.table (edges = key_nc, seconds = t_nc [["elapsed"]]), f_cache)
        message (sprintf ("network_cycles: %d cycles in %.1f s", length (key_nc), t_nc [["elapsed"]]))
    }

    # ---- wkpool::find_faces ------------------------------------------------------
    v <- unique (rbind (e [, .(id = node0, x = x0, y = y0)], e [, .(id = node1, x = x1, y = y1)]), by = "id")
    v [, .vx := seq_len (.N)]
    pool <- new_wkpool (
        vertices = data.frame (.vx = v$.vx, x = v$x, y = v$y),
        vx0 = v$.vx [match (e$node0, v$id)], vx1 = v$.vx [match (e$node1, v$id)]
    )
    t_ff <- system.time (ff <- find_faces (pool))
    ff <- as.data.table (ff)
    ff [, edge_ := e$edge_id [.segment]]
    # A face can have islands inside it (a loop on a dead-end road, a
    # roundabout): its walk goes round their outsides too, joined by
    # bridge edges. Split each face's non-bridge half-edges into rings
    # (connected pieces) and keep the outer ring, the one with positive
    # area; islands are walked clockwise (negative area). network_cycles()
    # returns outer rings only, so outer rings are what is compared.
    bounded <- ff [.face_area > 0 & !.bridge]
    xy <- data.table (.vx = v$.vx, x = v$x - mean (range (v$x)), y = v$y - mean (range (v$y)))
    bounded [, `:=` (ax = xy$x [.vx0], ay = xy$y [.vx0], bx = xy$x [.vx1], by = xy$y [.vx1])]
    gb <- igraph::graph_from_data_frame (bounded [, .(as.character (.vx0), as.character (.vx1))], directed = FALSE)
    memb <- igraph::components (gb)$membership
    bounded [, ring := memb [match (as.character (.vx0), names (memb))]]
    bounded [, ring_area := sum (ax * by - bx * ay) / 2, by = .(.face, ring)]
    n_rings <- bounded [, .(n = uniqueN (ring)), by = .face]
    outer <- bounded [ring_area > 0]
    key_ff <- outer [, .(key = paste (sort (unique (edge_)), collapse = ",")), by = .face]
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
        w <- outer [.face == fid]
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
    summary <- c (extra, list (
        edges = ne, vertices = nv,
        euler_V_minus_E_plus_F = nv - ne + nf, genus = (2 - (nv - ne + nf)) / 2,
        network_cycles = length (key_nc), network_cycles_duplicates = dup_nc,
        network_cycles_seconds = t_nc [["elapsed"]],
        find_faces_bounded = nrow (key_ff), find_faces_seconds = t_ff [["elapsed"]],
    find_faces_with_islands = sum (n_rings$n > 1),
        in_both = length (both), only_network_cycles = length (only_nc),
        only_network_cycles_outer = sum (only_nc %in% key_outer),
        only_network_cycles_union_of_faces = sum (only_nc_union),
        only_find_faces = length (only_ff),
        only_find_faces_self_crossing = sum (ff_only_crossing),
        find_faces_self_crossing_all = sum (all_crossing)
    ))
    jsonlite::write_json (summary, file.path (outdir, paste0 (prefix, "compare-summary.json")), auto_unbox = TRUE, pretty = TRUE)
    fwrite (data.table (method = "network_cycles", edges = only_nc), file.path (outdir, paste0 (prefix, "only-network-cycles.csv")))
    fwrite (data.table (method = "find_faces", edges = only_ff), file.path (outdir, paste0 (prefix, "only-find-faces.csv")))
    str (summary)
    invisible (summary)
}
