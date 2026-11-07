# city.infra

Working repo for applying open urban-analysis tooling to Hobart, Tasmania,
with an eye to contributions that have the most positive social impact in
collaboration with the [UrbanAnalyst](https://github.com/UrbanAnalyst) (UA)
project.

The questions driving this:

- Where can existing tools (hypertidy and related packages) close real gaps in
  the UA stack, by framing, by wiring, or by new development?
- What does it take to run a UA-style analysis on Hobart using authoritative
  local data (theLIST) alongside the global sources UA relies on?
- How do locally curated descriptions of infrastructure (accessibility,
  quality, lived experience) feed into network analysis?

See [PLAN.md](PLAN.md) for the current plan and status.

## Layout

| Path | What it holds |
|------|---------------|
| `PLAN.md` | The working plan, steps and status |
| `docs/` | Write-ups produced by each plan step |
| `index.qmd`, `mapper.qmd`, `_quarto.yml` | Quarto website, published to GitHub Pages |
| `site/mapper/` | Browser segment mapper (JS): `core.js` logic with tests, `app.js` UI |
| `apps/cycle_mapper/` | Original Shiny app for selecting and grouping theLIST network segments |
| `data/` | Provenance of data the apps and docs use |
| `notes/` | Parked material: links, ideas, meeting notes, loose ends |

## Related tools

- [wkpool](https://github.com/hypertidy/wkpool): vertex pool topology for
  wk geometry (segments, shared vertices, cycles, neighbours)
- [cdtr](https://github.com/hypertidy/cdtr): constrained Delaunay
  triangulation with refinement
- [traipse](https://github.com/AustralianAntarcticDivision/traipse): geodesic
  track metrics (distance, bearing, turn, speed, interpolation)
- [thelist-opendata](https://github.com/hypertidy/thelist-opendata):
  cloud-native mirror of LIST open data, Tasmania
- [splh.transport releases](https://github.com/mdsumner/splh.transport/releases/):
  Parquet cache of statewide Tasmania transport segments used by the app

## Website

The site is a Quarto project. Render locally with `quarto render` (output in
`_site/`); `.github/workflows/publish.yml` tests the mapper logic
(`node --test site/mapper/test/*.test.mjs`) and publishes to GitHub Pages on every push
to `main`.
