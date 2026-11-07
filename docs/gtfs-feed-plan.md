---
title: "Plan: a robust Tasmanian GTFS feed"
date: 2026-11-07
---

Status: **plan only, nothing set up yet.** Decisions for Michael are in the
last section.

## Why

The Hobart access run ([02-hobart-pipeline.md](02-hobart-pipeline.md),
friction point 5) had to use a GTFS copy from a third-party GitHub repository
(`devdennish/TransitTasmania`, feed version 20, valid 2026-06-03 to
2026-08-31) because the official sites were unreachable from the sandbox. That
copy has already expired, so the result cannot be re-run against the same
timetable from any official source, and there is no current feed either.

GTFS feeds are short-lived (a validity window of about three months, replaced
at each timetable change) and producers do not keep old versions online. Any
access result that is meant to be reproducible or compared over time needs its
own archive. The goal is:

1. **Get** the current official feed on a schedule, from more than one route.
2. **Version** it: an immutable, dated snapshot per distinct feed, never
   overwritten.
3. **Validate** it, and warn before the current feed expires.
4. **Archive** it next to the LIST mirror on source.coop, with a catalogue
   that answers "which snapshot covers date X?".
5. Make analysis scripts take a **snapshot id**, not a loose zip.

## What is out there (checked 2026-11-07)

| Source | URL | Status from here |
|---|---|---|
| Public Transport Tasmania statewide feed (official) | on transport.tas.gov.au; exact download URL to confirm | site unreachable from sandbox |
| Metro Tasmania, Hobart (official producer URL per Mobility Database) | `http://www.metrotas.com.au/wp-content/uploads/transit/Hobart/google_transit.zip` | site unreachable from sandbox |
| Metro Tasmania, Launceston / Burnie | same pattern, `.../transit/Launceston/...`, `.../transit/Burnie/...` | unreachable |
| SkyBus Hobart (official, separate feed) | Azure blob URL with a signed token valid to 2027-02-28 (in the Mobility Database catalogue) | not tried |
| Mobility Database "latest" mirror, source 665 | `https://storage.googleapis.com/storage/v1/b/mdb-latest/o/au-tasmania-metrotas-hobart-gtfs-665.zip?alt=media` | **reachable** |
| Mobility Database historical datasets | API at api.mobilitydatabase.org, needs a free refresh token | unreachable from sandbox |
| Transitland archive | transit.land, needs a free API key | unreachable from sandbox |
| Third-party GitHub copies | e.g. `devdennish/TransitTasmania` | last resort only |

What the Mobility Database mirror actually holds right now: the **statewide**
feed (agencies Metro Tasmania, Tassielink, Derwent Valley Link, Redline
Coaches, Derwent Ferries; 2,033 stops), calendar 2026-05-03 to 2026-07-25,
last fetched 2026-06-04. So its copy has also expired, and it has not picked up
the newer feed version 20. That strongly suggests the producer URL moved
around mid-2026, from the Metro Tasmania site to the Public Transport
Tasmania site, and the catalogue has not caught up. Two lessons:

- Do not rely on one URL. Keep a small list of source URLs, try them in order,
  and treat "nothing new for longer than the current feed's validity" as an
  alert, not as success.
- Reporting the new URL to the Mobility Database catalogue (a one-file PR to
  `MobilityData/mobility-database-catalogs`) is a cheap public good, and it
  makes their mirror a working fallback for everyone.

Small data-quality note already visible: `agency_timezone` is
`Australia/Sydney`, not `Australia/Hobart`. The two agree on offsets and DST
today, so it is harmless for routing, but a validator will report it and it is
worth passing back to the producer.

## Design

### 1. Fetch (scheduled)

A small script, run **weekly** (feeds change at timetable changes, not
monthly; weekly also gives four chances inside any expiry warning window):

- Read `sources.csv`: `source_id, priority, url, kind (official | mirror |
  third_party), notes`. Initial rows: Public Transport Tasmania statewide
  (official, URL to confirm), Metro Tas Hobart producer URL (official), MDB
  latest 665 (mirror), SkyBus (official, separate `source_id`).
- For each source in priority order, do a conditional GET using the stored
  `ETag` / `Last-Modified`. Stop at the first official source that answers.
- Check the zip is a real GTFS (has `agency.txt`, `stops.txt`, `trips.txt`,
  `stop_times.txt`, and `calendar.txt` or `calendar_dates.txt`).
- Compute a **content hash**: sha256 over the sorted list of
  `(file name, sha256 of file bytes)`. Zip timestamps change on every
  re-export, so the zip's own hash would create false "new" snapshots.
- If the content hash matches the latest snapshot, record the check in the
  fetch log and stop. Otherwise create a new snapshot.

### 2. Version (snapshot layout)

Interim home, as Michael suggested: a `gtfs/` prefix inside the existing
[tas-transport-web](https://source.coop/mdsumner/tas-transport-web) product,
following the same `latest.json` plus dated-prefix convention the segment
layer uses:

```
gtfs/
  catalog.csv                         one row per snapshot (see below)
  catalog.parquet                     same, for DuckDB / hyparquet
  <source_id>/
    latest.json                       {"snapshot": "2026-11-09", ...}
    <YYYY-MM-DD>/                     date first seen (UTC fetch date)
      gtfs.zip                        byte-for-byte as downloaded
      snapshot.json                   provenance + summary (below)
      validation/report.json          MobilityData validator output
      validation/report.html
      tables/                         optional derived Parquet (step 5)
```

- `source_id` is the feed, not the URL: `tas-statewide`, `skybus-hobart`.
  When an agency's URL changes, the `source_id` stays and `sources.csv`
  changes.
- Snapshot id is the date the content was **first seen**. A second distinct
  feed on the same day gets a suffix (`2026-11-09b`). Never overwrite a prefix
  once written.
- `snapshot.json` records: source url and kind, fetched_at, HTTP `ETag` and
  `Last-Modified`, zip sha256, content hash, `feed_info` fields if present
  (`feed_version`, `feed_start_date`, `feed_end_date`), service date range
  computed from `calendar` and `calendar_dates`, agency ids and names, row
  counts per table, validator error and warning counts, and the fetcher's git
  commit.
- `catalog.csv` columns: `source_id, snapshot, service_start, service_end,
  feed_version, content_sha256, source_kind, n_errors, n_warnings`. This is the
  index that answers "which snapshot covers 2026-08-04?".

source.coop has no object versioning, so immutability is by convention: the
writer refuses to upload into an existing prefix, and only `latest.json` and
the catalogue are ever rewritten.

### 3. Validate and alert

- Run the
  [MobilityData GTFS validator](https://github.com/MobilityData/gtfs-validator)
  (Java CLI, writes JSON and HTML). Store the report with the snapshot. A feed
  with errors is **still archived** (it is what was published) but flagged in
  the catalogue.
- Fast checks in R with `gtfsrouter::extract_gtfs()` and
  `gtfs_timetable()` on one weekday, so we know the feed actually routes with
  the tools the analysis uses.
- **Expiry warning:** if the newest snapshot's `service_end` is within 21 days,
  or already past, the workflow opens (or updates) a GitHub issue. This is the
  check that would have caught the mid-2026 URL move.
- Fetch failures on every official source also open an issue, with the HTTP
  status per URL.

### 4. Backfill

Seed the archive with what can still be found, each labelled by `source_kind`:

- Mobility Database latest 665 (statewide, 2026-05-03 to 2026-07-25), kind
  `mirror`.
- The `devdennish/TransitTasmania` copy used by the Hobart run (version 20,
  2026-06-03 to 2026-08-31), kind `third_party`, so the published result is
  reproducible from our own store.
- Older versions from the Mobility Database historical API or Transitland, if
  Michael is happy to register for a free key.

### 5. Derived tables (later, optional)

Per snapshot, each GTFS table as Parquet, and `stops` as GeoParquet in
EPSG:4326 (for the browser mapper) with a 28355 copy for joining to LIST. This
lets the segment mapper show stops and lets DuckDB query timetables across
snapshots without unzipping. Not needed for the first version.

### 6. Use in analysis

- `scripts/02-hobart/03_access_rhh.R` takes a snapshot id (or "latest")
  instead of a zip path, resolves it via `latest.json` / `catalog.csv`, and
  picks the analysis date automatically: the first Tuesday inside the service
  range that is not in `calendar_dates` as a removed-service exception.
- Every output records `gtfs_source_id` and `gtfs_snapshot`, the same way
  saved mapper groupings record the segment snapshot.
- Re-running against a sequence of snapshots gives access change over time
  for free, which is a strong equity story (service cuts and additions by
  suburb).

## Where it runs

The open question is reachability. The sandbox block was a proxy allowlist,
not necessarily a geo-block, so a GitHub Actions runner may reach the
official sites fine. Order of preference:

1. **GitHub Actions**, weekly cron plus manual dispatch, uploading with the
   same `SOURCE_COOP_*` service-account secrets planned for the
   thelist-opendata `refresh.yml`. A one-off probe job (curl each source URL,
   print status and headers, upload nothing) settles reachability.
2. If official sites block cloud runners: run the same script from a machine
   in Tasmania (cron on Michael's machine or a small AU-region VM), with the
   Actions job falling back to the mirror only.

Code location: start in this repo as `scripts/gtfs/` plus
`.github/workflows/gtfs-fetch.yml`, kept self-contained (one R or Python
script, `sources.csv`, no dependency on the rest of city.infra) so it can move
into thelist-opendata or a future org repo unchanged. Tooling: R with `httr2`,
`digest`, `gtfsrouter` and `arrow`, or plain Python with `requests`; plus the
validator jar. Uploads via `aws s3 cp` against `https://data.source.coop`.

## Licence

To confirm before publishing: the terms on the Public Transport Tasmania GTFS
page (Tasmanian Government open data is usually CC BY 4.0 or 3.0 AU) and the
SkyBus feed terms. Record them in `tas-transport-web`'s README and in each
`snapshot.json`. Third-party copies are archived for reproducibility only and
labelled as such.

## Proposed steps

1. Michael confirms the current official download URL(s) from a browser in
   Tasmania, and the licence wording.
2. Add the probe workflow; run it once to settle where fetching can run.
3. Add the fetch / validate / upload script and weekly workflow, writing to
   `tas-transport-web/gtfs/` (needs the source.coop secrets on this repo).
4. Backfill the two known 2026 copies; re-run the Hobart access script
   against `gtfs/tas-statewide/<snapshot>` and against the current feed.
5. PR the updated producer URL to the Mobility Database catalogue; tell the
   producer about `agency_timezone`.
6. Later: derived Parquet tables, stops in the mapper, GTFS-Realtime if
   Tasmania publishes it, and a permanent home once the data layout settles.

## Decisions for Michael

- **Source URL:** what is the current official statewide GTFS download link?
- **Runner:** OK to add a probe-only Actions job first (no uploads)?
- **Scope:** statewide feed only to start, or also SkyBus Hobart?
- **Secrets:** add `SOURCE_COOP_*` credentials with write access to
  `tas-transport-web` on this repo (or run uploads from thelist-opendata)?
- **Backfill:** register for a Mobility Database or Transitland key to pull
  older versions, or start from the two copies we have?
