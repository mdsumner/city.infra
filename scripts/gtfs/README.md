# GTFS snapshots (plan step 1)

Fetch Tasmanian GTFS feeds and keep immutable, dated snapshots. Design and
open decisions: [docs/gtfs-feed-plan.md](../../docs/gtfs-feed-plan.md).

This folder is self-contained (no dependency on the rest of city.infra) so it
can move to thelist-opendata or a future org repo unchanged.

| File | What |
|---|---|
| `sources.csv` | Feeds (`source_id`) and the URLs to try for each, in priority order |
| `gtfs_snapshot.R` | Fetch, check, content-hash, dedup, write the store |

Needs R with httr2, digest, jsonlite and data.table; gtfsrouter for `--check`;
Java 17 and the [MobilityData validator](https://github.com/MobilityData/gtfs-validator)
CLI jar for `--validator`.

## Run

```sh
# fetch every feed in sources.csv into a local store
Rscript gtfs_snapshot.R --store gtfs-store --check \
    --validator gtfs-validator-8.0.1-cli.jar [--java /path/to/java]

# add a zip you already have (backfill)
Rscript gtfs_snapshot.R --store gtfs-store --file old.zip \
    --source-id metrotas-hobart --kind third_party --url <where it came from> \
    --snapshot 2026-06-03 --fetched-at 2026-06-03
```

Per feed, URLs are tried in order. The first that gives a valid feed (or a
`304 Not Modified`) ends the search. A feed is a new snapshot only if its
**content hash** (sha256 over the sorted file names and per-file sha256)
is new for that `source_id`. Re-zipped copies of the same feed are
recognised as duplicates.

## Store layout

Same layout as the planned `tas-transport-web/gtfs/` prefix:

```
catalog.csv                     one row per snapshot: service range, hashes, source
status.json, status.md          per feed: newest service_end, days left, ok / expiring /
                                expired, validator counts, this run's fetch attempts,
                                and alert reasons
<source_id>/latest.json         snapshot with the latest service_end
<source_id>/fetch_log.csv       every check, including unchanged / not modified / failed
<source_id>/<YYYY-MM-DD>/gtfs.zip       byte-for-byte as downloaded
<source_id>/<YYYY-MM-DD>/snapshot.json  provenance, HTTP validators, agencies,
                                        per-file sha256 + row counts, routing check,
                                        validator summary
<source_id>/<YYYY-MM-DD>/validation/    validator report.json, report.html
```

The snapshot id is the UTC date the content was first seen (with a `b`, `c`
suffix for a second feed on the same day). The script refuses to write into an
existing snapshot directory.

An alert is raised (`status.json` `alert: true`) when a feed is expired or
within 21 days of expiry, when every URL for a feed fails, or when a feed came
from a mirror because the official URL failed (the sign that a producer URL
moved). Validator errors are reported but do not raise an alert.

## Workflows

- `.github/workflows/gtfs-weekly.yml`: Mondays early (Hobart time), plus manual
  runs. Fetches, validates and checks every feed, writes `status.md` to the
  run summary, keeps the store as a 90-day run artifact, and keeps one open
  "GTFS feed alert" issue in sync: created on the first alert, body updated
  each run, a comment when what is wrong changes, closed when all is OK.
- `.github/workflows/gtfs-probe.yml`: manual, probe-only reachability check.

Not done yet: upload to source.coop. Until then nothing persists between
weekly runs, so every run shows each feed as "new" and the run artifacts
(90 days) are the only copies.
