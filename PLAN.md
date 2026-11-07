# Plan

Goal: make this suite (wkpool, cdtr, traipse, the patch / cycle mapper app,
and the theLIST data mirror) a dependable, collaborator-ready resource for
UrbanAnalyst-style work, grounded in Hobart and its local data, with the
best positive social impact for the effort.

## Done so far

- **Project home.** This repo, restructured into `apps/`, `docs/`, `notes/`,
  with a Quarto site and a browser segment mapper (`mapper.qmd`) on LISTmap
  basemap tiles.
- **UA stack map.** What each active UA repo does and where our tools plug in:
  [docs/01-ua-stack.md](docs/01-ua-stack.md). Cycle mapper review:
  [docs/cycle-mapper-review.md](docs/cycle-mapper-review.md).
- **Data.** theLIST mirror on source.coop (`mdsumner/thelist-opendata`) and a
  derived web layer (`mdsumner/tas-transport-web`, EPSG:4326 Parquet +
  PMTiles). Status and loose ends:
  [notes/thelist-opendata-status.md](notes/thelist-opendata-status.md).
- **First Hobart run.** Walk + public transport time to the Royal Hobart
  Hospital for every Greater Hobart address, with dodgr / gtfsrouter on the
  LIST network, plus a friction log:
  [docs/02-hobart-pipeline.md](docs/02-hobart-pipeline.md).
- **sfnetworks review.** Node matching by rounded-coordinate strings is behind
  sfnetworks issue #130; wkpool already matches on exact coordinate values
  (PR #7, `notes/sfnetworks/sfnetworks-wkpool.md`).
- **GTFS feed.** Planning a dependable source of dated Tasmanian GTFS
  snapshots (interim home: `tas-transport-web`, `gtfs/` prefix).

## Steps

1. **Make the data dependable** - not started.
   Merge the thelist-opendata refresh workflow
   (hypertidy/thelist-opendata#1, needs `SOURCE_COOP_*` secrets), add a
   README (licence, provenance, snapshot layout) to `tas-transport-web`, and
   land dated GTFS snapshots. Outcome: anyone can reproduce the Hobart runs
   from public, versioned inputs.

2. **Turn friction into upstream-ready fixes** - not started.
   Draft here first, nothing filed upstream without sign-off:
   - a LIST road-class (`TRAN_CLASS`) weighting profile for dodgr;
   - non-OSM network input and a configurable travel-time cap (currently
     fixed at 60 min) for m4ra;
   - the wkpool-based node matching fix for sfnetworks #130, with an R
     benchmark.
   Output: `docs/03-upstream-candidates.md` plus draft issues / patches.

3. **Write the suite's collaborator intro** - not started.
   One Quarto page per tool (wkpool, cdtr, traipse, patch app): what it does,
   shown on Hobart data, and which UA problem it addresses.
   Output: `docs/suite/` pages on the site.

4. **Shortlist 2 to 3 joint proposals** - not started.
   Scored by social impact against effort, including where the patch app
   sits alongside UA `neighbourhoods` and `superblock`. Each with a concrete
   first step.
   Output: `docs/04-proposals.md`

Steps 1 and 2 can run in parallel.

## Conventions

- Publishing uses Quarto (website on GitHub Pages).
- Do not name collaborators in docs or PRs for now.

## Open decisions

- Name for a future GitHub org for this work (ideas so far: `netpatch`,
  `placegrain`, `urbanpatch`; availability not checked).
- Long-term home for GTFS snapshots and the derived web layer.
