# Walk and walk + public transport time from every Greater Hobart address
# to the Royal Hobart Hospital (48 Liverpool St), weekday morning.
#
# Inputs (paths as arguments, see README.md in this folder):
#   edges.csv      from 01_extract_network.py (LIST transport segments)
#   addresses.csv  from 02_extract_addresses.py (LIST address points)
#   gtfs.zip       Tasmanian GTFS feed (Public Transport Tasmania)
#   outdir         where results are written
#
# Usage: Rscript 03_access_rhh.R edges.csv addresses.csv gtfs.zip outdir

suppressPackageStartupMessages({
    library (dodgr)
    library (gtfsrouter)
    library (m4ra)
    library (data.table)
    library (sf)
})

args <- commandArgs (trailingOnly = TRUE)
if (length (args) != 4L) {
    stop ("Usage: Rscript 03_access_rhh.R edges.csv addresses.csv gtfs.zip outdir")
}
f_edges <- args [1]
f_addr <- args [2]
f_gtfs <- args [3]
outdir <- args [4]
dir.create (outdir, showWarnings = FALSE, recursive = TRUE)

walk_speed <- 4.8 / 3.6 # m/s
rhh_geocode <- 278074L # LIST GEOCODE_ID for 48 Liverpool St, Hobart
gtfs_date <- 20260804L # a Tuesday inside the feed's validity (2026-06-03 to 2026-08-31)
gtfs_day <- "tuesday"
start_window <- c (7.5, 8.5) * 3600 # first boarding between 07:30 and 08:30
max_pt_seconds <- 2 * 3600 # cap on the public transport leg
use_m4ra <- identical (Sys.getenv ("USE_M4RA"), "1") # see friction log

tick <- function (msg) message (format (Sys.time (), "%H:%M:%S"), "  ", msg)

# ---- 1. walk graph from LIST segments ------------------------------------

tick ("reading LIST segments")
e <- fread (f_edges)
not_walkable <- c ("Railway", "Railway Siding", "Ferry", "Tramway")
e <- e [TRANS_TYPE != "Rail" & !TRAN_CLASS %in% not_walkable]
e <- e [node0 != node1] # self-loops carry no routing information

one_way <- data.frame (
    edge_id = paste0 (e$TRANSEG_ID, "f"),
    from_id = as.character (e$node0),
    from_lon = e$x0, from_lat = e$y0,
    to_id = as.character (e$node1),
    to_lon = e$x1, to_lat = e$y1,
    d = e$COMP_LEN,
    transeg_id = e$TRANSEG_ID,
    tran_class = e$TRAN_CLASS,
    lga = sub (";.*", "", e$lga)
)
rev_way <- one_way
rev_way$edge_id <- paste0 (e$TRANSEG_ID, "r")
rev_way [, c ("from_id", "from_lon", "from_lat", "to_id", "to_lon", "to_lat")] <-
    one_way [, c ("to_id", "to_lon", "to_lat", "from_id", "from_lon", "from_lat")]
graph <- rbind (one_way, rev_way)
graph$d_weighted <- graph$d # walking: no preference by class yet (see friction log)
graph <- dodgr_components (graph)
n_comp <- length (unique (graph$component))
graph <- graph [graph$component == 1L, ]
tick (sprintf (
    "walk graph: %d directed edges, %d vertices (largest of %d components)",
    nrow (graph), nrow (dodgr_vertices (graph)), n_comp
))
v <- dodgr_vertices (graph)
v_lga <- unique (rbind (
    data.frame (id = graph$from_id, lga = graph$lga),
    data.frame (id = graph$to_id, lga = graph$lga)
))
v_lga <- v_lga [!duplicated (v_lga$id), ]

# ---- 2. address points ---------------------------------------------------

tick ("reading LIST address points")
a <- fread (f_addr)
a_sf <- st_as_sf (a, coords = c ("EASTING", "NORTHING"), crs = 28355)
a_ll <- st_coordinates (st_transform (a_sf, 4326))
a$lon <- a_ll [, 1]
a$lat <- a_ll [, 2]
bb <- c (range (v$x), range (v$y))
a <- a [lon >= bb [1] & lon <= bb [2] & lat >= bb [3] & lat <= bb [4]]
a$vid <- v$id [match_points_to_verts (v, a [, c ("lon", "lat")])]
vi <- match (a$vid, v$id)
a$snap_m <- geodist::geodist (
    a [, c ("lon", "lat")], v [vi, c ("x", "y")],
    paired = TRUE, measure = "haversine"
)
a <- a [snap_m <= 250]
a$lga <- v_lga$lga [match (a$vid, v_lga$id)]
tick (sprintf ("%d addresses within 250 m of the walk graph", nrow (a)))

rhh <- a [GEOCODE_ID == rhh_geocode]
stopifnot (nrow (rhh) == 1L)

# ---- 3. walk only --------------------------------------------------------

tick ("walk-only times")
to_v <- unique (a$vid)
d_walk <- dodgr_dists (graph, from = rhh$vid, to = to_v) [1, ]
a$walk_min <- d_walk [match (a$vid, to_v)] / walk_speed / 60

# ---- 4. public transport, stop to stop -----------------------------------

tick ("reading GTFS")
gtfs <- extract_gtfs (f_gtfs, quiet = TRUE)
stops_in <- gtfs$stops [
    stop_lon >= bb [1] & stop_lon <= bb [2] & stop_lat >= bb [3] & stop_lat <= bb [4]
]
gtfs <- gtfs_transfer_table (gtfs, d_limit = 200, quiet = TRUE)
gtfs <- gtfs_timetable (gtfs, day = gtfs_day, date = gtfs_date)
tick (sprintf (
    "%d stops in the study area of %d in the feed; computing stop-to-stop matrix",
    nrow (stops_in), nrow (gtfs$stops)
))
if (use_m4ra) {
    # m4ra's stop matrix: fast, but max_traveltime is fixed at 60 minutes
    tt <- m4ra_gtfs_traveltimes (
        gtfs,
        start_time_limits = start_window,
        day = gtfs_day,
        from_stops = stops_in$stop_id,
        next_interval = FALSE,
        quiet = TRUE
    )
    dur <- tt$duration # seconds; rows = from_stops, cols = all stops
} else {
    # Same matrix from gtfsrouter's public API, with a longer time cap
    hms_to_s <- function (x) {
        p <- strsplit (x, ":", fixed = TRUE)
        vapply (p, function (i) sum (as.numeric (i) * c (3600, 60, 1)), numeric (1))
    }
    rows <- parallel::mclapply (stops_in$stop_id, function (s) {
        r <- gtfs_traveltimes (gtfs, from = s, start_time_limits = start_window,
                               from_is_id = TRUE, max_traveltime = max_pt_seconds,
                               quiet = TRUE)
        out <- rep (NA_real_, nrow (gtfs$stops))
        if (nrow (r) > 0) out [match (r$stop_id, gtfs$stops$stop_id)] <- hms_to_s (r$duration)
        out
    }, mc.cores = parallel::detectCores ())
    dur <- do.call (rbind, rows)
    dimnames (dur) <- list (stops_in$stop_id, gtfs$stops$stop_id)
}
tick (sprintf ("stop matrix %d x %d", nrow (dur), ncol (dur)))

# ---- 5. stop -> hospital, then address -> hospital ------------------------

stop_v <- v$id [match_points_to_verts (v, gtfs$stops [, c ("stop_lon", "stop_lat")])]
walk_stop_rhh <- dodgr_dists (graph, from = rhh$vid, to = unique (stop_v)) [1, ]
walk_stop_rhh <- walk_stop_rhh [match (stop_v, unique (stop_v))] / walk_speed # s, all stops

col_stop <- match (colnames (dur), gtfs$stops$stop_id)
row_stop <- match (rownames (dur), gtfs$stops$stop_id)
egress <- walk_stop_rhh [col_stop]
via_pt <- apply (dur, 1, function (r) suppressWarnings (min (r + egress, na.rm = TRUE)))
via_pt [!is.finite (via_pt)] <- NA
# best time from boarding stop s to the hospital door
t_stop <- pmin (via_pt, walk_stop_rhh [row_stop], na.rm = TRUE)

# Multi-source trick: one virtual source joined to each boarding stop's
# vertex by an edge whose length is the remaining time to the hospital
# expressed as walking distance. Walk graph is symmetric, so the
# shortest path source -> address equals address -> stop -> hospital.
ok <- !is.na (t_stop)
src <- data.frame (
    edge_id = paste0 ("src", seq_len (sum (ok))),
    from_id = "SRC", from_lon = rhh$lon, from_lat = rhh$lat,
    to_id = stop_v [row_stop [ok]],
    to_lon = v$x [match (stop_v [row_stop [ok]], v$id)],
    to_lat = v$y [match (stop_v [row_stop [ok]], v$id)],
    d = t_stop [ok] * walk_speed,
    transeg_id = NA, tran_class = "virtual", lga = NA, component = 1L
)
src$d_weighted <- src$d
src <- src [, names (graph)]
g2 <- rbind (graph, src)
tick ("walk + public transport times")
d_pt <- dodgr_dists (g2, from = "SRC", to = to_v) [1, ]
a$pt_min <- pmin (d_pt [match (a$vid, to_v)] / walk_speed / 60, a$walk_min, na.rm = TRUE)

# ---- 6. summaries ---------------------------------------------------------

band <- function (x) cut (x, c (0, 30, 45, 60, Inf), right = TRUE,
                          labels = c ("<=30", "30-45", "45-60", ">60"))
s <- a [, .(
    addresses = .N,
    walk_le30 = round (100 * mean (walk_min <= 30, na.rm = TRUE), 1),
    pt_le30 = round (100 * mean (pt_min <= 30, na.rm = TRUE), 1),
    pt_le45 = round (100 * mean (pt_min <= 45, na.rm = TRUE), 1),
    pt_le60 = round (100 * mean (pt_min <= 60, na.rm = TRUE), 1),
    pt_median_min = round (median (pt_min, na.rm = TRUE))
), by = lga] [order (pt_median_min)]
tot <- a [, .(
    lga = "all", addresses = .N,
    walk_le30 = round (100 * mean (walk_min <= 30, na.rm = TRUE), 1),
    pt_le30 = round (100 * mean (pt_min <= 30, na.rm = TRUE), 1),
    pt_le45 = round (100 * mean (pt_min <= 45, na.rm = TRUE), 1),
    pt_le60 = round (100 * mean (pt_min <= 60, na.rm = TRUE), 1),
    pt_median_min = round (median (pt_min, na.rm = TRUE))
)]
s <- rbind (s, tot)
fwrite (s, file.path (outdir, "rhh-access-by-lga.csv"))
print (s)

loc <- a [, .(addresses = .N, pt_median_min = round (median (pt_min, na.rm = TRUE)),
              lga = lga [1]), by = LOCALITY] [addresses >= 200] [order (-pt_median_min)]
fwrite (loc, file.path (outdir, "rhh-access-by-locality.csv"))
print (head (loc, 15))

fwrite (a [, .(GEOCODE_ID, LOCALITY, lga, lon, lat, snap_m,
               walk_min = round (walk_min, 1), pt_min = round (pt_min, 1))],
        file.path (outdir, "rhh-access-addresses.csv.gz"))

png (file.path (outdir, "rhh-access-map.png"), width = 1600, height = 1600, res = 200)
cols <- hcl.colors (4, "viridis", rev = TRUE)
b <- band (a$pt_min)
plot (a$lon, a$lat, pch = ".", col = cols [as.integer (b)], asp = 1 / cos (-42.88 * pi / 180),
      xlab = "", ylab = "", axes = FALSE,
      main = "Walk + public transport to Royal Hobart Hospital\nweekday, first boarding 07:30-08:30 (minutes)")
points (rhh$lon, rhh$lat, pch = 4, cex = 2, lwd = 2, col = "red")
legend ("topright", legend = levels (b), col = cols, pch = 15, bty = "n")
mtext ("Address points and transport segments: theLIST, State of Tasmania. GTFS: Public Transport Tasmania",
       side = 1, cex = 0.6)
invisible (dev.off ())
tick ("done")
