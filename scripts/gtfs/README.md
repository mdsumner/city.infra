# GTFS snapshots (plan step 1)

Fetch Tasmanian GTFS feeds and keep immutable, dated snapshots. Design and
open decisions: [docs/gtfs-feed-plan.md](../../docs/gtfs-feed-plan.md).

This folder is self-contained (no dependency on the rest of city.infra) so it
can move to thelist-opendata or a future org repo unchanged.

| File | What |
|---|---|
| `sources.csv` | Feeds (`source_id`) and the URLs to try for each, in priority order |
| `gtfs_snapshot.R` | Fetch, check, content-hash, dedup, write the store |

Needs R with httr2, digest, jsonlite and data.table; gtfsrouter for `--check`.

## Run

```sh
# fetch every feed in sources.csv into a local store
Rscript gtfs_snapshot.R --store gtfs-store --check

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
status.json                     newest service_end per feed, days left, ok / expiring / expired
<source_id>/latest.json         snapshot with the latest service_end
<source_id>/fetch_log.csv       every check, including unchanged / not modified / failed
<source_id>/<YYYY-MM-DD>/gtfs.zip       byte-for-byte as downloaded
<source_id>/<YYYY-MM-DD>/snapshot.json  provenance, HTTP validators, agencies,
                                        per-file sha256 + row counts, routing check
```

The snapshot id is the UTC date the content was first seen (with a `b`, `c`
suffix for a second feed on the same day). The script refuses to write into an
existing snapshot directory.

Not done yet: MobilityData validator, upload to source.coop, the weekly
workflow and the expiry issue. `.github/workflows/gtfs-probe.yml` is a
manual, probe-only job that checks which URLs a GitHub runner can reach.
