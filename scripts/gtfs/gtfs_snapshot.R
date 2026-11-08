#!/usr/bin/env Rscript
# Fetch Tasmanian GTFS feeds and keep immutable, dated snapshots.
# Design: docs/gtfs-feed-plan.md. Self-contained: needs only httr2, digest,
# jsonlite, data.table (and gtfsrouter for --check).
#
# Fetch every feed in sources.csv (weekly job):
#   Rscript gtfs_snapshot.R --store DIR [--sources sources.csv] [--only ID] [--check]
#
# Add a zip you already have (backfill):
#   Rscript gtfs_snapshot.R --store DIR --file x.zip --source-id ID \
#       --kind third_party --url URL [--snapshot YYYY-MM-DD] [--check]
#
# Store layout (mirrors tas-transport-web/gtfs/):
#   DIR/catalog.csv                      one row per snapshot
#   DIR/status.json                      newest service_end per feed, days left
#   DIR/<source_id>/latest.json
#   DIR/<source_id>/fetch_log.csv        every check, including "unchanged"
#   DIR/<source_id>/<YYYY-MM-DD>/gtfs.zip, snapshot.json
#
# A snapshot directory is never written twice. Only latest.json, catalog.csv,
# status.json and fetch_log.csv are rewritten.

suppressPackageStartupMessages ({
    library (data.table)
})

warn_days <- 21L
required <- c ("agency.txt", "stops.txt", "routes.txt", "trips.txt", "stop_times.txt")

# ---- args -------------------------------------------------------------------

parse_args <- function (a) {
    out <- list ()
    i <- 1L
    while (i <= length (a)) {
        k <- sub ("^--", "", a [i])
        if (k %in% c ("check")) {
            out [[k]] <- TRUE
            i <- i + 1L
        } else {
            out [[gsub ("-", "_", k)]] <- a [i + 1L]
            i <- i + 2L
        }
    }
    out
}

script_dir <- function () {
    f <- sub ("^--file=", "", grep ("^--file=", commandArgs (FALSE), value = TRUE))
    if (length (f)) dirname (normalizePath (f)) else getwd ()
}

# ---- reading a feed ---------------------------------------------------------

sha256_file <- function (f) digest::digest (file = f, algo = "sha256")

gtfs_date <- function (x) as.Date (as.character (x), format = "%Y%m%d")

#' Summarise a GTFS zip. Errors if it is not a usable GTFS.
summarise_gtfs <- function (zip) {
    td <- tempfile ("gtfs")
    dir.create (td)
    on.exit (unlink (td, recursive = TRUE))
    files <- utils::unzip (zip, exdir = td, junkpaths = TRUE)
    nm <- basename (files)
    missing <- setdiff (required, nm)
    if (!any (c ("calendar.txt", "calendar_dates.txt") %in% nm)) {
        missing <- c (missing, "calendar.txt|calendar_dates.txt")
    }
    if (length (missing)) stop ("not a GTFS feed, missing: ", paste (missing, collapse = ", "))

    # Content hash: independent of zip timestamps and compression.
    o <- order (nm)
    file_sha <- vapply (files [o], sha256_file, "")
    content_sha <- digest::digest (paste0 (nm [o], " ", file_sha, "\n", collapse = ""),
        algo = "sha256", serialize = FALSE)

    rd <- function (f) fread (file.path (td, f), colClasses = "character", encoding = "UTF-8")
    rows <- vapply (nm [o], function (f) nrow (rd (f)), 0L)

    starts <- ends <- as.Date (character ())
    if ("calendar.txt" %in% nm) {
        cal <- rd ("calendar.txt")
        starts <- c (starts, gtfs_date (cal$start_date))
        ends <- c (ends, gtfs_date (cal$end_date))
    }
    if ("calendar_dates.txt" %in% nm) {
        cd <- rd ("calendar_dates.txt")
        added <- gtfs_date (cd$date [cd$exception_type == "1"])
        starts <- c (starts, added)
        ends <- c (ends, added)
    }
    feed_info <- if ("feed_info.txt" %in% nm) as.list (rd ("feed_info.txt") [1]) else NULL
    ag <- rd ("agency.txt")

    list (
        content_sha256 = content_sha,
        service_start = format (min (starts)),
        service_end = format (max (ends)),
        feed_info = feed_info,
        agencies = ag [, intersect (c ("agency_id", "agency_name", "agency_timezone"), names (ag)), with = FALSE],
        files = data.table (file = nm [o], sha256 = unname (file_sha), rows = unname (rows))
    )
}

#' Can gtfsrouter build a timetable for a weekday inside the service range?
routing_check <- function (zip, s) {
    if (!requireNamespace ("gtfsrouter", quietly = TRUE)) {
        return (list (ok = NA, note = "gtfsrouter not installed"))
    }
    d0 <- max (as.Date (s$service_start), Sys.Date ())
    if (d0 > as.Date (s$service_end)) d0 <- as.Date (s$service_start)
    days <- seq (d0, by = "day", length.out = 7)
    d <- days [format (days, "%u") == "2"] # Tuesday
    res <- tryCatch ({
        g <- gtfsrouter::extract_gtfs (zip, quiet = TRUE)
        g <- gtfsrouter::gtfs_timetable (g, day = "tuesday", date = as.integer (format (d, "%Y%m%d")))
        list (ok = nrow (g$timetable) > 0, date = format (d), timetable_rows = nrow (g$timetable))
    }, error = function (e) list (ok = FALSE, date = format (d), error = conditionMessage (e)))
    res
}

# ---- store ------------------------------------------------------------------

catalog_cols <- c ("source_id", "snapshot", "service_start", "service_end", "feed_version",
    "content_sha256", "source_kind", "source_url", "fetched_at", "n_errors", "n_warnings")

read_catalog <- function (store) {
    f <- file.path (store, "catalog.csv")
    if (file.exists (f)) {
        fread (f, colClasses = "character", na.strings = "")
    } else {
        setnames (data.table (matrix (character (), 0, length (catalog_cols))), catalog_cols)
    }
}

new_snapshot_id <- function (src_dir, date) {
    id <- format (date)
    for (sfx in c ("", letters [-1])) {
        cand <- paste0 (id, sfx)
        if (!dir.exists (file.path (src_dir, cand))) return (cand)
    }
    stop ("too many snapshots on ", id)
}

git_commit <- function () {
    r <- tryCatch (system2 ("git", c ("-C", script_dir (), "rev-parse", "HEAD"),
        stdout = TRUE, stderr = FALSE), error = function (e) NA_character_, warning = function (w) NA_character_)
    if (length (r)) r [1] else NA_character_
}

write_json <- function (x, f) {
    jsonlite::write_json (x, f, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null", digits = NA)
}

#' Add a zip to the store, unless its content is already there.
#' Returns the snapshot id (new or existing).
add_snapshot <- function (store, zip, source_id, kind, url, fetched_at = Sys.time (),
                          http = list (), snapshot = NULL, check = FALSE) {
    s <- summarise_gtfs (zip)
    cat_dt <- read_catalog (store)
    sid <- source_id
    dup <- cat_dt [cat_dt$source_id == sid & cat_dt$content_sha256 == s$content_sha256]
    if (nrow (dup)) {
        message (source_id, ": content already archived as ", dup$snapshot [1])
        return (invisible (list (snapshot = dup$snapshot [1], new = FALSE)))
    }

    src_dir <- file.path (store, source_id)
    dir.create (src_dir, recursive = TRUE, showWarnings = FALSE)
    id <- if (is.null (snapshot)) new_snapshot_id (src_dir, as.Date (fetched_at, tz = "UTC")) else snapshot
    snap_dir <- file.path (src_dir, id)
    if (dir.exists (snap_dir)) stop ("refusing to write into existing snapshot ", snap_dir)
    dir.create (snap_dir)
    file.copy (zip, file.path (snap_dir, "gtfs.zip"))

    meta <- list (
        source_id = source_id,
        snapshot = id,
        source_kind = kind,
        source_url = url,
        fetched_at = format (fetched_at, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
        http = http,
        zip_sha256 = sha256_file (zip),
        zip_bytes = file.size (zip),
        content_sha256 = s$content_sha256,
        service_start = s$service_start,
        service_end = s$service_end,
        feed_info = s$feed_info,
        agencies = s$agencies,
        files = s$files,
        validation = NULL, # MobilityData validator: not wired up yet
        routing_check = if (check) routing_check (zip, s) else NULL,
        fetcher_commit = git_commit ()
    )
    write_json (meta, file.path (snap_dir, "snapshot.json"))

    row <- data.table (
        source_id = source_id, snapshot = id,
        service_start = s$service_start, service_end = s$service_end,
        feed_version = if (is.null (s$feed_info$feed_version)) NA_character_ else s$feed_info$feed_version,
        content_sha256 = s$content_sha256, source_kind = kind, source_url = url,
        fetched_at = meta$fetched_at, n_errors = NA_character_, n_warnings = NA_character_
    )
    cat_dt <- rbind (cat_dt, row, fill = TRUE)
    setorder (cat_dt, source_id, snapshot)
    fwrite (cat_dt, file.path (store, "catalog.csv"))
    update_latest (store, source_id, cat_dt)
    message (source_id, ": new snapshot ", id, " (service ", s$service_start, " to ", s$service_end, ")")
    invisible (list (snapshot = id, new = TRUE))
}

#' latest.json points at the snapshot whose service runs latest, so that a
#' backfilled older feed never displaces the current one.
update_latest <- function (store, source_id, cat_dt = read_catalog (store)) {
    want <- source_id
    x <- cat_dt [cat_dt$source_id == want]
    x <- x [order (service_end, snapshot, decreasing = TRUE)] [1]
    write_json (list (
        source_id = source_id, snapshot = x$snapshot,
        service_start = x$service_start, service_end = x$service_end,
        path = paste0 (source_id, "/", x$snapshot, "/gtfs.zip")
    ), file.path (store, source_id, "latest.json"))
}

# ---- fetching ---------------------------------------------------------------

log_fetch <- function (store, source_id, row) {
    f <- file.path (store, source_id, "fetch_log.csv")
    dir.create (dirname (f), recursive = TRUE, showWarnings = FALSE)
    fwrite (as.data.table (row), f, append = file.exists (f))
}

last_validators <- function (store, source_id, url) {
    f <- file.path (store, source_id, "fetch_log.csv")
    if (!file.exists (f)) return (list ())
    lg <- fread (f, colClasses = "character")
    want <- url
    lg <- lg [lg$url == want & status %in% c ("200", "304") & outcome != "invalid"]
    if (!nrow (lg)) return (list ())
    as.list (lg [.N, .(etag, last_modified)])
}

# ETags are logged without their quotes (CSV-safe) and re-quoted when sent.
etag_strip <- function (e) gsub ('"', "", e, fixed = TRUE)
etag_wrap <- function (e) {
    if (startsWith (e, "W/")) paste0 ('W/"', substring (e, 3), '"') else paste0 ('"', e, '"')
}

fetch_url <- function (url, validators = list ()) {
    dest <- tempfile (fileext = ".zip")
    req <- httr2::request (url) |>
        httr2::req_user_agent ("city.infra-gtfs-archiver (+https://github.com/mdsumner/city.infra)") |>
        httr2::req_timeout (120) |>
        httr2::req_retry (max_tries = 3) |>
        httr2::req_error (is_error = function (r) FALSE)
    if (length (validators$etag) && !is.na (validators$etag) && nzchar (validators$etag)) {
        req <- httr2::req_headers (req, `If-None-Match` = etag_wrap (validators$etag))
    }
    if (length (validators$last_modified) && !is.na (validators$last_modified) && nzchar (validators$last_modified)) {
        req <- httr2::req_headers (req, `If-Modified-Since` = validators$last_modified)
    }
    resp <- tryCatch (httr2::req_perform (req, path = dest), error = function (e) e)
    if (inherits (resp, "error")) {
        return (list (status = NA_integer_, error = conditionMessage (resp)))
    }
    hdr <- function (h) {
        v <- httr2::resp_header (resp, h)
        if (is.null (v)) NA_character_ else v
    }
    list (
        status = httr2::resp_status (resp), path = dest, final_url = resp$url,
        etag = etag_strip (hdr ("etag")), last_modified = hdr ("last-modified"), error = NA_character_
    )
}

# Never write signed-URL secrets (Azure SAS "sig", AWS "Signature") to logs.
redact_url <- function (u) {
    if (is.null (u) || is.na (u)) return (u)
    gsub ("((^|[?&])(sig|Signature|X-Amz-Signature)=)[^&]*", "\\1REDACTED", u)
}

#' A source row with resolve = "mdb_catalog" gives the raw URL of a Mobility
#' Database catalogue entry; the download URL is read from it at fetch time,
#' so renewed signed URLs (e.g. SkyBus) are picked up without editing here.
resolve_url <- function (url, resolve) {
    if (is.na (resolve) || !nzchar (resolve)) return (url)
    if (resolve != "mdb_catalog") stop ("unknown resolve: ", resolve)
    j <- jsonlite::fromJSON (url)
    j$urls$direct_download
}

#' Try each URL for one feed in priority order. Stops at the first URL that
#' gives a valid feed or says "not modified".
fetch_source <- function (store, src, check = FALSE) {
    src <- src [order (as.integer (priority))]
    sid <- src$source_id [1]
    for (i in seq_len (nrow (src))) {
        u <- src$url [i]
        res_kind <- if ("resolve" %in% names (src)) src$resolve [i] else NA_character_
        r <- tryCatch (fetch_url (resolve_url (u, res_kind), last_validators (store, sid, u)),
            error = function (e) list (status = NA_integer_, error = conditionMessage (e)))
        r$final_url <- redact_url (r$final_url)
        row <- list (checked_at = format (Sys.time (), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
            url = u, kind = src$kind [i], status = as.character (r$status),
            etag = if (is.null (r$etag)) NA else r$etag,
            last_modified = if (is.null (r$last_modified)) NA else r$last_modified,
            outcome = NA_character_, snapshot = NA_character_, error = r$error)
        if (isTRUE (r$status == 304L)) {
            row$outcome <- "not_modified"
        } else if (isTRUE (r$status == 200L)) {
            res <- tryCatch (
                add_snapshot (store, r$path, sid, src$kind [i], u,
                    http = list (final_url = r$final_url, etag = r$etag, last_modified = r$last_modified),
                    check = check),
                error = function (e) e)
            if (inherits (res, "error")) {
                row$outcome <- "invalid"
                row$error <- conditionMessage (res)
            } else {
                row$outcome <- if (res$new) "new" else "unchanged"
                row$snapshot <- res$snapshot
            }
        } else {
            row$outcome <- "failed"
        }
        log_fetch (store, sid, row)
        message (sid, " <- ", u, ": ", row$status, " ", row$outcome)
        if (row$outcome %in% c ("new", "unchanged", "not_modified")) return (invisible (row))
    }
    invisible (NULL)
}

#' Per-feed summary for alerting: newest service_end and days left.
write_status <- function (store) {
    cat_dt <- read_catalog (store)
    if (!nrow (cat_dt)) return (invisible ())
    today <- Sys.Date ()
    st <- cat_dt [, .SD [order (service_end, snapshot, decreasing = TRUE)] [1], by = source_id]
    st [, days_left := as.integer (as.Date (service_end) - today)]
    st [, state := fifelse (days_left < 0, "expired", fifelse (days_left <= warn_days, "expiring", "ok"))]
    out <- st [, .(source_id, snapshot, service_start, service_end, days_left, state, source_kind)]
    write_json (list (checked = format (today), warn_days = warn_days, feeds = out),
        file.path (store, "status.json"))
    print (out)
    invisible (out)
}

# ---- main -------------------------------------------------------------------

main <- function (a = parse_args (commandArgs (TRUE))) {
    if (is.null (a$store)) stop ("--store DIR is required")
    dir.create (a$store, recursive = TRUE, showWarnings = FALSE)
    check <- isTRUE (a$check)

    if (!is.null (a$file)) {
        for (k in c ("source_id", "kind", "url")) {
            if (is.null (a [[k]])) stop ("--file needs --", gsub ("_", "-", k))
        }
        add_snapshot (a$store, a$file, a$source_id, a$kind, a$url,
            fetched_at = if (is.null (a$fetched_at)) Sys.time () else as.POSIXct (a$fetched_at, tz = "UTC"),
            snapshot = a$snapshot, check = check)
    } else {
        f <- if (is.null (a$sources)) file.path (script_dir (), "sources.csv") else a$sources
        srcs <- fread (f, colClasses = "character")
        if (!is.null (a$only)) srcs <- srcs [source_id %in% strsplit (a$only, ",") [[1]]]
        for (sid in unique (srcs$source_id)) {
            if (is.null (fetch_source (a$store, srcs [source_id == sid], check))) {
                message (sid, ": every source failed")
            }
        }
    }
    write_status (a$store)
}

if (sys.nframe () == 0L) main ()
