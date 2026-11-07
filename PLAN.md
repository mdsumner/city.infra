# Plan

Goal: identify and start on the contributions to UrbanAnalyst-style analysis
that give the best positive social impact, grounded in Hobart and its local
data, and using the tools we already have where they fit.

## Steps

1. **Map the UA stack** - in progress.
   What each active UA repo does, what it takes in and gives out, where the
   R / C++ / Rust boundaries sit, and where wkpool, cdtr, traipse, theLIST
   data and the cycle mapper app could plug in.
   Output: [docs/01-ua-stack.md](docs/01-ua-stack.md)

1b. **Cycle mapper to a static Quarto site** - review done, port not started.
   Review and approach: [docs/cycle-mapper-review.md](docs/cycle-mapper-review.md).
   Network cache moves into thelist-opendata on S3; the mapper becomes a
   Quarto OJS page on GitHub Pages alongside the other resources.

2. **Run Hobart through it** - not started.
   Inventory local data (theLIST transport and cadastre, OSM, Metro Tasmania
   GTFS, ABS census and SEIFA, population grids, elevation) and attempt a
   minimal end-to-end pipeline (dodgr / m4ra to an access result), logging
   every point of friction.
   Output: `docs/02-hobart-pipeline.md` plus scripts.

3. **Match tools to gaps** - not started.
   Position the cycle mapper alongside UA `neighbourhoods` and `superblock`,
   and sort all candidate contributions into framing, wiring and new
   development, scored by social impact against effort.
   Output: `docs/03-candidates.md`

4. **Shortlist 2 to 3 proposals** - not started.
   Each with a concrete first step, ready to take to the UA maintainers.
   Output: `docs/04-proposals.md`

## Conventions

- Publishing uses Quarto (website on GitHub Pages).

## Open decisions

- Name for a future GitHub org for this work (ideas so far: `netpatch`,
  `placegrain`, `urbanpatch`; availability not checked).
