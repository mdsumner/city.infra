# Cycle Infrastructure Segment Mapper

Map cycle infrastructure descriptions to theLIST Tasmania transport network.

## Quick start

```r
install.packages(c("shiny", "leaflet", "dplyr", "wk", "arrow", "DT", "PROJ"))
shiny::runApp("cycle_mapper")
```

## Search syntax

- `Argyle` - finds all segments with "Argyle" in name
- `Argyle [Federal, Burnett]` - Argyle between Federal and Burnett intersections
- `Macquarie [Darcy, Elboden], Weld, Anglesea [Adelaide, Macquarie]` - multi-segment route

Comma separates route parts. Each part can be a whole street or a bounded section.

## Workflow

1. **Search** - enter street name (with optional cross-street bounds)
2. **Select** - click linestrings on map (green = selected)
3. **Save** - give it an ID and description (auto-saves to `cycle_mappings.csv`)
4. **Progress tab** - see all saved mappings on a map

## Output format

`cycle_mappings.csv`:
```
id,description,transeg_id
CP001,"Argyle St bike lane",12345
CP001,"Argyle St bike lane",12346
```

## Features

- Segments ordered spatially (walk the network)
- Cross-street filtering to find specific sections
- Multi-segment route definitions
- Context segments shown (3 either side of bounded section)
- Auto-save/load - persists between sessions
- Progress map showing all saved mappings
