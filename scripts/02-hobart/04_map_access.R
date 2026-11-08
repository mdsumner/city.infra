# Map walk + public transport time to the Royal Hobart Hospital, with
# geographic context: land / water, LGA boundaries, suburb names, the River
# Derwent and kunanyi / Mt Wellington.
#
# Inputs:
#   addresses   rhh-access-addresses.csv.gz from 03_access_rhh.R
#   context     folder with theLIST GeoParquet (see README.md):
#                 lga.parquet          list_local_govt_areas (statewide)
#                 loc_<lga>.parquet    list_locality_postcode, one per LGA
#   run.json    from 03_access_rhh.R (analysis date, feed)
#   out.png
#
# Usage: Rscript 04_map_access.R addresses.csv.gz context/ run.json out.png

suppressPackageStartupMessages ({
    library (data.table)
    library (sf)
})

args <- commandArgs (trailingOnly = TRUE)
if (length (args) != 4L) {
    stop ("Usage: Rscript 04_map_access.R addresses.csv.gz context/ run.json out.png")
}
f_addr <- args [1]
ctx <- args [2]
run <- jsonlite::fromJSON (args [3])
f_out <- args [4]

crs <- 28355 # GDA94 / MGA zone 55, metres

read_gpq <- function (f) {
    x <- as.data.frame (arrow::read_parquet (f))
    g <- sf::st_as_sfc (structure (x$geometry, class = "WKB"), crs = crs)
    x$geometry <- NULL
    x$geometry_bbox <- NULL
    sf::st_sf (x, geometry = g)
}

# ---- data ------------------------------------------------------------------

a <- fread (f_addr)
a <- a [!is.na (pt_min)]
xy <- st_coordinates (st_transform (st_as_sf (a, coords = c ("lon", "lat"), crs = 4326), crs))
a [, `:=` (x = xy [, 1], y = xy [, 2])]

ll <- function (lon, lat) st_coordinates (st_transform (st_sfc (st_point (c (lon, lat)), crs = 4326), crs))

# Main map: the urban area, Brighton to Margate and Richmond to Lauderdale.
# The overview inset shows the whole study area.
p0 <- ll (147.17, -43.05)
p1 <- ll (147.55, -42.63)
xlim <- c (p0 [1], p1 [1])
ylim <- c (p0 [2], p1 [2])
pad <- 2500
xlim_all <- range (a$x, xlim) + c (-pad, pad)
ylim_all <- range (a$y, ylim) + c (-pad, pad)
frame <- st_as_sfc (st_bbox (c (xmin = min (xlim, xlim_all) - 20000, ymin = min (ylim, ylim_all) - 20000,
    xmax = max (xlim, xlim_all) + 20000, ymax = max (ylim, ylim_all) + 20000), crs = crs))

lga <- read_gpq (file.path (ctx, "lga.parquet"))
lga <- suppressWarnings (st_intersection (st_make_valid (lga), frame))
loc <- do.call (rbind, lapply (list.files (ctx, "^loc_.*parquet$", full.names = TRUE), function (f) {
    x <- read_gpq (f)
    x [, c ("NAME", "geometry")]
}))
# Per-LGA files can each hold a clipped piece of the same locality
loc$NAME <- toupper (loc$NAME)
loc <- aggregate (loc ["NAME"], by = list (name = loc$NAME), FUN = function (x) x [1])
loc$NAME <- loc$name

# ---- styling ---------------------------------------------------------------

# Ordinal blue ramp, validated (dataviz validate_palette.js --ordinal): light
# = quick, dark = slow. Land is near-white so the light step still reads;
# water is a darker neutral, never blue, so it is not mistaken for data.
bands <- c (0, 30, 45, 60, Inf)
band_lab <- c ("up to 30", "30 to 45", "45 to 60", "over 60")
ramp <- c ("#86b6ef", "#3987e5", "#1c5cab", "#0d366b")
col_land <- "#fcfcfb"
col_water <- "#dedcd6"
col_lga <- "#8a8984"
ink <- "#1f1f1d"
ink_muted <- "#5f5e5a"

a [, band := cut (pt_min, bands, right = TRUE, labels = band_lab)]
a <- a [order (band)] # slow on top, so outlying slow areas are not hidden

# Suburbs to name: the larger ones plus places discussed in the write-up
named <- c ("SANDY BAY", "NEW TOWN", "MOONAH", "GLENORCHY", "CLAREMONT",
    "AUSTINS FERRY", "GRANTON", "BRIDGEWATER", "GAGEBROOK", "BRIGHTON", "OLD BEACH",
    "LINDISFARNE", "BELLERIVE", "HOWRAH", "ROKEBY", "LAUDERDALE", "CAMBRIDGE",
    "RICHMOND", "SANDFORD", "SOUTH ARM", "KINGSTON", "BLACKMANS BAY", "HUNTINGFIELD",
    "MARGATE", "SNUG", "KETTERING", "TAROONA", "LENAH VALLEY", "RISDON VALE",
    "CLARENDON VALE", "SORELL", "MIDWAY POINT", "OPOSSUM BAY", "FERN TREE")
lab <- loc [loc$NAME %in% named, ]
lab_xy <- st_coordinates (st_point_on_surface (st_geometry (lab)))
in_frame <- lab_xy [, 1] > xlim [1] & lab_xy [, 1] < xlim [2] &
    lab_xy [, 2] > ylim [1] & lab_xy [, 2] < ylim [2]
lab <- lab [in_frame, ]
lab_xy <- lab_xy [in_frame, , drop = FALSE]
title_case <- function (x) gsub ("\\b([A-Z])([A-Z]*)", "\\1\\L\\2", x, perl = TRUE)

rhh <- a [GEOCODE_ID == 278074]
if (!nrow (rhh)) rhh <- data.table (x = ll (147.3297, -42.8813) [1], y = ll (147.3297, -42.8813) [2])
kunanyi <- ll (147.2366, -42.8961) # summit

# ---- draw ------------------------------------------------------------------

asp_h <- diff (ylim) / diff (xlim)
w_map <- 10   # main map width, inches
w_side <- 3.2 # side panel: legend, overview, scale
h_map <- w_map * asp_h
png (f_out, width = w_map + w_side, height = h_map + 1.2, units = "in", res = 300)
par (bg = col_land, family = "sans", oma = c (0.8, 0, 3.2, 0))
layout (matrix (1:2, 1), widths = c (w_map, w_side))

# ---- 1. main map ------------------------------------------------------------
par (mar = c (0.3, 0.3, 0.3, 0.3))
plot.new ()
plot.window (xlim, ylim, asp = 1, xaxs = "i", yaxs = "i")
usr <- par ("usr")
rect (usr [1], usr [3], usr [2], usr [4], col = col_water, border = NA)
plot (st_geometry (lga), col = col_land, border = NA, add = TRUE)
points (a$x, a$y, pch = 15, cex = 0.12, col = ramp [as.integer (a$band)])
plot (st_geometry (lga), col = NA, border = col_lga, lwd = 0.6, lty = "22", add = TRUE)

# labels with a halo so they read over the points
halo_text <- function (x, y, labels, col = ink, cex = 0.62, font = 1, halo = col_land, r = 0.6, adj = NULL) {
    d <- r * strwidth ("x", cex = cex) * 0.25
    for (dx in c (-d, 0, d)) for (dy in c (-d, 0, d)) {
        if (dx != 0 || dy != 0) text (x + dx, y + dy, labels, col = halo, cex = cex, font = font, adj = adj)
    }
    text (x, y, labels, col = col, cex = cex, font = font, adj = adj)
}
halo_text (lab_xy [, 1], lab_xy [, 2], title_case (lab$NAME))

rv <- ll (147.36, -42.95)
halo_text (rv [1], rv [2], "River Derwent", col = ink_muted, cex = 0.85, font = 3, halo = col_water)
points (kunanyi [1], kunanyi [2], pch = 17, cex = 1.1, col = ink)
halo_text (kunanyi [1], kunanyi [2] - 1100, "kunanyi /\nMt Wellington", cex = 0.62, font = 2)
# hospital: marker, leader line to a label on the water
lx <- rhh$x + 3200
ly <- rhh$y - 2600
segments (rhh$x, rhh$y, lx - 150, ly + 300, col = ink, lwd = 0.8)
points (rhh$x, rhh$y, pch = 21, cex = 1.6, lwd = 2, bg = "#ffffff", col = "#c0362c")
halo_text (lx, ly, "Royal Hobart Hospital\n(destination)", cex = 0.72, font = 2, halo = col_water, adj = c (0, 1))
box (col = col_lga, lwd = 0.6)
upi <- diff (usr [1:2]) / par ("pin") [1] # map metres per inch, for the scale bar

# ---- 2. side panel ------------------------------------------------------------
par (mar = c (0.3, 0.5, 0.3, 0.3))
plot.new ()
pin <- par ("pin")
plot.window (c (0, pin [1] * upi), c (0, pin [2] * upi), xaxs = "i", yaxs = "i") # same metres per inch
u <- par ("usr")
yy <- function (f) u [3] + f * diff (u [3:4]) # vertical position as a fraction

# legend
n <- a [, .N, by = band] [order (band)]
pct <- round (100 * n$N [match (band_lab, n$band)] / sum (n$N))
text (u [1], yy (0.985), "Walk + public transport\nto the hospital door", adj = c (0, 1), cex = 0.85, font = 2, col = ink)
ky <- yy (0.935) - (0:3) * 0.028 * diff (u [3:4])
sw <- 0.07 * diff (u [1:2])
rect (u [1], ky - 0.009 * diff (u [3:4]), u [1] + sw, ky + 0.009 * diff (u [3:4]), col = ramp, border = NA)
text (u [1] + 1.4 * sw, ky, sprintf ("%s min", band_lab), adj = c (0, 0.5), cex = 0.75, col = ink)
text (u [2] - 0.02 * diff (u [1:2]), ky, sprintf ("%d%%", pct), adj = c (1, 0.5), cex = 0.75, col = ink_muted)
text (u [1], ky [4] - 0.03 * diff (u [3:4]), "Share of all addresses in the study area.\nFirst boarding 07:30 to 08:30;\nno waiting time counted.",
    adj = c (0, 1), cex = 0.6, col = ink_muted)

# overview of the whole study area
ov_top <- yy (0.70)
ov_w <- diff (u [1:2]) * 0.98
ov_h <- ov_w * diff (ylim_all) / diff (xlim_all)
ov_h <- min (ov_h, ov_top - yy (0.20))
ov_w <- ov_h * diff (xlim_all) / diff (ylim_all)
text (u [1], ov_top + 0.012 * diff (u [3:4]), "Whole study area (box = main map)", adj = c (0, 0), cex = 0.7, font = 2, col = ink)
sx_ <- function (x) u [1] + (x - xlim_all [1]) / diff (xlim_all) * ov_w
sy_ <- function (y) ov_top - ov_h + (y - ylim_all [1]) / diff (ylim_all) * ov_h
rect (u [1], ov_top - ov_h, u [1] + ov_w, ov_top, col = col_water, border = NA)
lga_ov <- suppressWarnings (st_intersection (lga, st_as_sfc (st_bbox (c (xmin = xlim_all [1], ymin = ylim_all [1],
    xmax = xlim_all [2], ymax = ylim_all [2]), crs = crs))))
for (g in st_geometry (lga_ov)) {
    for (pg in if (inherits (g, "MULTIPOLYGON")) g else list (g)) {
        polygon (sx_ (pg [[1]] [, 1]), sy_ (pg [[1]] [, 2]), col = col_land, border = col_lga, lwd = 0.3)
    }
}
points (sx_ (a$x), sy_ (a$y), pch = 15, cex = 0.03, col = ramp [as.integer (a$band)])
rect (sx_ (xlim [1]), sy_ (ylim [1]), sx_ (xlim [2]), sy_ (ylim [2]), border = ink, lwd = 1)
points (sx_ (rhh$x), sy_ (rhh$y), pch = 21, cex = 0.6, lwd = 1, bg = "#ffffff", col = "#c0362c")
rect (u [1], ov_top - ov_h, u [1] + ov_w, ov_top, border = col_lga, lwd = 0.6)

# scale bar at the main map's scale, and north arrow
sx <- u [1] + 400
sy <- yy (0.10)
rect (sx, sy, sx + 2500, sy + 250, col = ink, border = NA)
rect (sx + 2500, sy, sx + 5000, sy + 250, col = col_land, border = ink, lwd = 0.8)
text (c (sx, sx + 2500, sx + 5000), sy - 500, c ("0", "2.5", "5 km"), cex = 0.62, col = ink, adj = c (0.5, 1))
text (u [1], sy + 800, "Main map scale", adj = c (0, 0), cex = 0.62, col = ink_muted)
nx <- u [2] - 0.12 * diff (u [1:2])
arrows (nx, sy + 1500, nx, sy + 4000, length = 0.08, lwd = 1.5, col = ink)
text (nx, sy + 4700, "N", cex = 0.8, font = 2, col = ink)

# ---- titles -------------------------------------------------------------------
d <- as.Date (as.character (run$gtfs_date), "%Y%m%d")
mtext ("Time from home to the Royal Hobart Hospital by walking and public transport",
    side = 3, outer = TRUE, line = 1.6, adj = 0.005, cex = 1.05, font = 2, col = ink)
mtext (sprintf ("%s. Each dot is an address (%s). GTFS feed version %s.",
    paste (format (d, "%A"), as.integer (format (d, "%d")), format (d, "%B %Y")), format (nrow (a), big.mark = ","),
    if (is.null (run$gtfs_feed_version)) "?" else run$gtfs_feed_version),
    side = 3, outer = TRUE, line = 0.4, adj = 0.005, cex = 0.72, col = ink_muted)
mtext ("Address points, transport segments, LGA and locality boundaries: theLIST, State of Tasmania (CC BY). GTFS: Public Transport Tasmania (CC BY 4.0). Projection: GDA94 / MGA zone 55.",
    side = 1, outer = TRUE, line = 0, adj = 0.005, cex = 0.5, col = ink_muted)
invisible (dev.off ())
message ("wrote ", f_out)
