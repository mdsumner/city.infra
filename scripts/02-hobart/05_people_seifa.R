# People and disadvantage: weight the address-level access results by 2021
# Census population (mesh blocks) and attach SEIFA 2021 IRSD by SA1.
#
# Inputs (ABS, see README.md for URLs):
#   addresses   rhh-access-addresses.csv.gz from 03_access_rhh.R
#   mb.zip      MB_2021_AUST_SHP_GDA94.zip (mesh block boundaries, ASGS Ed. 3)
#   counts      Mesh Block Counts, 2021.xlsx (Table 6 = Tasmania)
#   seifa       Statistical Area Level 1, Indexes, SEIFA 2021.xlsx (Table 2 = IRSD)
#   outdir
#
# Method: each address is placed in its 2021 mesh block; the mesh block's
# usual resident population is shared equally among the study-area
# addresses in it. Each address gets its SA1's IRSD decile within Tasmania
# (1 = most disadvantaged) and nationally.
#
# Usage: Rscript 05_people_seifa.R addresses.csv.gz mb.zip counts.xlsx seifa.xlsx outdir

suppressPackageStartupMessages ({
    library (data.table)
    library (sf)
    library (readxl)
})
sf_use_s2 (FALSE)

args <- commandArgs (trailingOnly = TRUE)
if (length (args) != 5L) stop ("Usage: Rscript 05_people_seifa.R addresses.csv.gz mb.zip counts.xlsx seifa.xlsx outdir")
f_addr <- args [1]
f_mb <- args [2]
f_counts <- args [3]
f_seifa <- args [4]
outdir <- args [5]
dir.create (outdir, showWarnings = FALSE, recursive = TRUE)
tick <- function (msg) message (format (Sys.time (), "%H:%M:%S"), "  ", msg)

# ---- inputs -------------------------------------------------------------------

tick ("addresses")
a <- fread (f_addr)
a <- a [!is.na (pt_min)]
a_sf <- st_as_sf (a, coords = c ("lon", "lat"), crs = 4283, remove = FALSE) # WGS84 ~ GDA94 at metre level

tick ("mesh blocks (Tasmania, intersecting the study area)")
bb <- st_bbox (a_sf)
wkt <- st_as_text (st_as_sfc (bb))
mb <- st_read (paste0 ("/vsizip/", f_mb, "/MB_2021_AUST_GDA94.shp"), quiet = TRUE,
    query = "SELECT MB_CODE21, MB_CAT21, SA1_CODE21, SA2_NAME21, GCC_CODE21 FROM MB_2021_AUST_GDA94 WHERE STE_CODE21 = '6'",
    wkt_filter = wkt)
mb <- st_make_valid (mb)

tick ("mesh block counts")
cnt <- as.data.table (read_excel (f_counts, sheet = "Table 6", skip = 6, .name_repair = "minimal"))
setnames (cnt, c ("MB_CODE21", "MB_CAT", "AREA_SQKM", "dwellings", "persons", "state") [seq_len (ncol (cnt))])
# ABS codes may come through as numbers or text, depending on the sheet
as_code <- function (x) if (is.numeric (x)) format (x, scientific = FALSE, trim = TRUE) else trimws (x)
cnt [, MB_CODE21 := as_code (MB_CODE21)]
cnt <- cnt [grepl ("^[0-9]{11}$", MB_CODE21)]
stopifnot (nrow (cnt) > 10000, all (as.integer (cnt$state) == 6L)) # Tasmania
cnt [, `:=` (persons = as.numeric (persons), dwellings = as.numeric (dwellings))]

tick ("SEIFA IRSD by SA1")
irsd <- as.data.table (read_excel (f_seifa, sheet = "Table 2", skip = 6, col_names = FALSE, .name_repair = "minimal"))
irsd <- irsd [, c (1, 2, 3, 6, 11)]
irsd [[1]] <- as_code (irsd [[1]])
setnames (irsd, c ("SA1_CODE21", "sa1_pop", "irsd_score", "irsd_decile_aus", "irsd_decile_tas"))
irsd <- irsd [grepl ("^[0-9]{11}$", SA1_CODE21)]
stopifnot (nrow (irsd) > 50000)
irsd [, `:=` (irsd_score = as.numeric (irsd_score), irsd_decile_aus = as.integer (irsd_decile_aus),
    irsd_decile_tas = as.integer (irsd_decile_tas))]

# ---- join ---------------------------------------------------------------------

tick ("addresses -> mesh blocks")
j <- st_join (a_sf, mb, join = st_within)
j <- as.data.table (st_drop_geometry (j))
j <- j [!duplicated (GEOCODE_ID)]
j <- merge (j, cnt [, .(MB_CODE21, persons, dwellings)], by = "MB_CODE21", all.x = TRUE)
j [, n_addr_mb := .N, by = MB_CODE21]
j [, people := fifelse (is.na (persons), 0, persons / n_addr_mb)]
j <- merge (j, irsd [, .(SA1_CODE21, irsd_score, irsd_decile_tas, irsd_decile_aus)], by = "SA1_CODE21", all.x = TRUE)

# Coverage: population of the mesh blocks that hold study addresses, against
# all residents of Greater Hobart mesh blocks
mb_hit <- unique (j$MB_CODE21)
gh <- cnt [MB_CODE21 %in% mb$MB_CODE21 [mb$GCC_CODE21 == "6GHOB"]]
coverage <- list (
    addresses = nrow (j),
    addresses_without_mb = sum (is.na (j$MB_CODE21)),
    people_allocated = round (sum (j$people)),
    greater_hobart_residents = sum (gh$persons),
    greater_hobart_residents_in_mbs_with_study_addresses = sum (gh$persons [gh$MB_CODE21 %in% mb_hit]),
    people_without_irsd = round (sum (j$people [is.na (j$irsd_decile_tas)]))
)
stopifnot (coverage$people_allocated > 0)
tick (sprintf ("%s people allocated; Greater Hobart residents %s, of whom %s live in mesh blocks with study addresses",
    format (coverage$people_allocated, big.mark = ","), format (coverage$greater_hobart_residents, big.mark = ","),
    format (coverage$greater_hobart_residents_in_mbs_with_study_addresses, big.mark = ",")))

# ---- summaries ----------------------------------------------------------------

wmedian <- function (x, w) {
    o <- order (x)
    x <- x [o]
    w <- w [o]
    x [which (cumsum (w) >= sum (w) / 2) [1]]
}
summ <- function (d) d [, .(
    people = round (sum (people)),
    addresses = .N,
    pt_le30 = round (100 * sum (people [pt_min <= 30]) / sum (people), 1),
    pt_le45 = round (100 * sum (people [pt_min <= 45]) / sum (people), 1),
    pt_le60 = round (100 * sum (people [pt_min <= 60]) / sum (people), 1),
    pt_median_min = round (wmedian (pt_min, people))
)]

j [, irsd_quintile_tas := (irsd_decile_tas + 1L) %/% 2L]
q_lab <- c ("1 most disadvantaged", "2", "3", "4", "5 least disadvantaged")

by_q <- j [people > 0 & !is.na (irsd_quintile_tas), summ (.SD), by = irsd_quintile_tas] [order (irsd_quintile_tas)]
by_q [, quintile := q_lab [irsd_quintile_tas]]
by_lga <- rbind (j [people > 0, summ (.SD), by = lga] [order (pt_median_min)],
    cbind (lga = "all", summ (j [people > 0])))
by_lga_q <- j [people > 0 & !is.na (irsd_quintile_tas) & lga %in% c ("hobart", "glenorchy", "clarence", "kingborough", "brighton"),
    summ (.SD), by = .(lga, irsd_quintile_tas)] [order (lga, irsd_quintile_tas)]
# Most disadvantaged areas with the longest trips: SA2s with >= 1000 people
by_sa2 <- j [people > 0 & !is.na (irsd_score), c (summ (.SD), .(
    irsd_decile_tas_median = as.integer (round (wmedian (irsd_decile_tas, people))),
    irsd_score_mean = round (weighted.mean (irsd_score, people)))),
    by = SA2_NAME21] [people >= 1000] [order (-pt_median_min)]

fwrite (by_q [, c ("quintile", setdiff (names (by_q), c ("quintile", "irsd_quintile_tas"))), with = FALSE],
    file.path (outdir, "rhh-people-by-irsd-quintile.csv"))
fwrite (by_lga, file.path (outdir, "rhh-people-by-lga.csv"))
fwrite (by_lga_q, file.path (outdir, "rhh-people-by-lga-irsd-quintile.csv"))
fwrite (by_sa2, file.path (outdir, "rhh-people-by-sa2.csv"))
jsonlite::write_json (coverage, file.path (outdir, "people-coverage.json"), auto_unbox = TRUE, pretty = TRUE)

# ---- chart: SA2 median time against disadvantage ------------------------------
# One series: each SA2 (>= 1000 people) is a dot sized by population. Lower
# IRSD score = more disadvantaged, on the left. Selective direct labels.
ink <- "#1f1f1d"
ink_muted <- "#5f5e5a"
grid <- "#e6e5e1"
surface <- "#fcfcfb"
dot <- "#2a78d6"
png (file.path (outdir, "rhh-people-sa2.png"), width = 10, height = 6.5, units = "in", res = 200)
par (mar = c (4.6, 4.6, 4.2, 1), bg = surface, family = "sans", col.axis = ink_muted, fg = ink_muted)
x <- by_sa2$irsd_score_mean
y <- by_sa2$pt_median_min
ylim <- c (0, max (y) * 1.05)
yt <- seq (0, ylim [2], by = 25)
plot (x, y, type = "n", xlab = "", ylab = "", ylim = ylim, axes = FALSE, yaxs = "i")
abline (h = yt, col = grid, lwd = 1)
axis (1, lwd = 0, lwd.ticks = 0.6, cex.axis = 0.8)
axis (2, at = yt, lwd = 0, las = 1, cex.axis = 0.8)
overall <- wmedian (j$pt_min [j$people > 0], j$people [j$people > 0])
abline (h = overall, col = ink_muted, lty = 2, lwd = 1)
text (par ("usr") [1], overall, sprintf ("median over all people: %d min", round (overall)), adj = c (0, -0.5), cex = 0.7, col = ink_muted)
cex <- 0.6 + 2.6 * sqrt (by_sa2$people / max (by_sa2$people))
o <- order (-by_sa2$people)
points (x [o], y [o], pch = 21, cex = cex [o], bg = dot, col = surface, lwd = 1.5)
show <- by_sa2$SA2_NAME21 %in% c ("Bridgewater - Gagebrook", "Brighton - Pontville", "Rokeby", "Risdon Vale",
    "Clarendon Vale", "South Arm", "Bruny Island - Kettering", "Margate - Snug", "Hobart", "Sandy Bay",
    "Kingston Beach - Blackmans Bay", "Glenorchy", "Cambridge", "Berriedale - Chigwell", "Sorell - Richmond") |
    (by_sa2$irsd_decile_tas_median <= 2)
# labels go left of dots in the right third (room), and Risdon Vale goes left of Rokeby
lab_pos <- ifelse (x > quantile (par ("usr") [1:2], 0.5) + 0.3 * diff (par ("usr") [1:2]), 2, 4)
lab_pos [by_sa2$SA2_NAME21 == "Risdon Vale"] <- 2
text (x [show], y [show], by_sa2$SA2_NAME21 [show], pos = lab_pos [show], offset = 0.7, cex = 0.62, col = ink)
mtext ("IRSD score, SEIFA 2021 (people-weighted mean of SA1s)", side = 1, line = 2.4, cex = 0.8, col = ink)
mtext ("<- more disadvantaged", side = 1, line = 3.4, adj = 0, cex = 0.7, col = ink_muted)
mtext ("less disadvantaged ->", side = 1, line = 3.4, adj = 1, cex = 0.7, col = ink_muted)
mtext ("Median minutes to the hospital", side = 2, line = 3.2, cex = 0.8, col = ink)
mtext ("Walk + public transport to the Royal Hobart Hospital, by SA2", side = 3, line = 2.4, adj = 0, cex = 1, font = 2, col = ink)
mtext ("Each dot is an SA2 with at least 1,000 residents, sized by 2021 Census population. Median is over people.",
    side = 3, line = 1.2, adj = 0, cex = 0.72, col = ink_muted)
invisible (dev.off ())

print (by_q)
print (by_lga)
print (head (by_sa2, 15))
tick ("done")
