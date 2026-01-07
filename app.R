library(shiny)
library(leaflet)
library(dplyr)
library(wk)

# --- Configuration ---
NETWORK_URL <- "https://github.com/mdsumner/splh.transport/releases/download/latest/LIST_TRANSPORT_SEGMENTS.parquet"
HOBART_CENTER <- c(-42.88, 147.33)
MAPPINGS_FILE <- "cycle_mappings.csv"

# Visual style
STYLE <- list(
  candidate = list(color = "#888888", weight = 3, opacity = 0.6),
  context = list(color = "#cccccc", weight = 2, opacity = 0.4),  # lighter for context
  selected = list(color = "#00cc44", weight = 6, opacity = 0.9),
  saved = list(color = "#3366cc", weight = 4, opacity = 0.8),
  anchor = list(color = "#ff7700", radius = 10)
)

#' Load saved mappings from CSV
load_mappings <- function(file = MAPPINGS_FILE) {
  if (!file.exists(file)) return(list())
  
  df <- read.csv(file, stringsAsFactors = FALSE)
  if (nrow(df) == 0) return(list())
  
  # Convert to list grouped by id
  split(df, df$id) |>
    lapply(function(x) {
      list(
        id = x$id[1],
        description = x$description[1],
        transeg_ids = as.character(x$transeg_id)
      )
    })
}

#' Save mappings to CSV (overwrites duplicate IDs)
save_mappings <- function(mappings, file = MAPPINGS_FILE) {
  if (length(mappings) == 0) {
    write.csv(data.frame(id = character(), description = character(), 
                         transeg_id = character()), file, row.names = FALSE)
    return()
  }
  
  out <- do.call(rbind, lapply(mappings, function(x) {
    data.frame(id = x$id, description = x$description, 
               transeg_id = x$transeg_ids, stringsAsFactors = FALSE)
  }))
  write.csv(out, file, row.names = FALSE)
}

# --- Data loading ---
load_network <- function(url = NETWORK_URL) {
  message("Loading transport network...")
  segments <- arrow::read_parquet(url)
  
  # Convert WKB to wk, check CRS by coordinate range
  segments$geom <- wk::wkb(segments$SHAPE)
  
  # Get centroids for quick spatial filtering
  coords <- wk::wk_coords(segments$geom)
  centroids <- coords |>
    group_by(feature_id) |>
    summarise(cx = mean(x), cy = mean(y), .groups = "drop")
  
  segments$cx <- centroids$cx
  segments$cy <- centroids$cy
  
  # Check if already lon/lat
  x_range <- range(segments$cx, na.rm = TRUE)
  if (x_range[1] > 100 && x_range[2] < 160) {
    message("Data appears to be in lon/lat")
    segments$clng <- segments$cx
    segments$clat <- segments$cy
  } else {
    message("Data needs reprojection - using PROJ")
    ll <- PROJ::proj_trans(cbind(segments$cx, segments$cy), 
                           source = "EPSG:7855", target = "EPSG:4326")
    segments$clng <- ll[, 1]
    segments$clat <- ll[, 2]
  }
  
  message(sprintf("Loaded %d segments", nrow(segments)))
  segments
}

# Convert a wkb geometry to lon/lat coordinate list for leaflet
geom_to_lonlat <- function(geom_wkb, from_crs = "EPSG:7855") {
  coords <- wk::wk_coords(geom_wkb)
  
  # Check if reprojection needed
  # lon/lat: x in 100-180, y in -50 to -10 for Australia
  # MGA Zone 55: x ~ 200000-800000, y ~ 5000000-6000000
  if (max(coords$x) < 180 && min(coords$x) > 100) {
    # Already lon/lat
    list(lng = coords$x, lat = coords$y)
  } else {
    ll <- PROJ::proj_trans(cbind(coords$x, coords$y), 
                           source = from_crs, target = "EPSG:4326")
    list(lng = ll[, 1], lat = ll[, 2])
  }
}

#' Order segments by network connectivity
#' Returns row indices in spatial sequence order
order_segments_by_network <- function(geoms, tol = 1) {
  n <- length(geoms)
  if (n <= 1) return(seq_len(n))
  

  # Get start/end points of each linestring
  endpoints <- lapply(seq_len(n), function(i) {
    coords <- wk::wk_coords(geoms[i])
    list(
      start = c(coords$x[1], coords$y[1]),
      end = c(coords$x[nrow(coords)], coords$y[nrow(coords)])
    )
  })
  
  # Build adjacency: which segments share an endpoint?
  # Two segments are adjacent if any endpoint is within tolerance
  adjacency <- vector("list", n)
  for (i in seq_len(n)) {
    adjacency[[i]] <- integer()
    for (j in seq_len(n)) {
      if (i == j) next
      # Check all 4 endpoint combinations
      pts_i <- c(endpoints[[i]]$start, endpoints[[i]]$end)
      pts_j <- c(endpoints[[j]]$start, endpoints[[j]]$end)
      
      for (pi in list(endpoints[[i]]$start, endpoints[[i]]$end)) {
        for (pj in list(endpoints[[j]]$start, endpoints[[j]]$end)) {
          dist <- sqrt((pi[1] - pj[1])^2 + (pi[2] - pj[2])^2)
          if (dist <= tol) {
            adjacency[[i]] <- c(adjacency[[i]], j)
            break
          }
        }
        if (j %in% adjacency[[i]]) break
      }
    }
  }
  
  # Find "end" segments (those with only 1 neighbor) - good starting points
  neighbor_counts <- sapply(adjacency, length)
  ends <- which(neighbor_counts == 1)
  
  # If no clear ends (loop?), just start from segment with fewest neighbors
  start_idx <- if (length(ends) > 0) ends[1] else which.min(neighbor_counts)[1]
  
  # Walk the network greedily
  visited <- logical(n)
  order <- integer(n)
  current <- start_idx
  
  for (step in seq_len(n)) {
    order[step] <- current
    visited[current] <- TRUE
    
    # Find unvisited neighbors
    neighbors <- adjacency[[current]]
    unvisited <- neighbors[!visited[neighbors]]
    
    if (length(unvisited) == 0) {
      # No connected unvisited - find nearest unvisited segment
      remaining <- which(!visited)
      if (length(remaining) == 0) break
      
      # Pick nearest by centroid distance
      current_centroid <- (endpoints[[current]]$start + endpoints[[current]]$end) / 2
      dists <- sapply(remaining, function(r) {
        r_centroid <- (endpoints[[r]]$start + endpoints[[r]]$end) / 2
        sqrt(sum((current_centroid - r_centroid)^2))
      })
      current <- remaining[which.min(dists)]
    } else {
      current <- unvisited[1]
    }
  }
  
  order[order > 0]  # trim if we broke early
}

#' Parse search syntax: 
#'   "MainStreet [Cross1, Cross2]" - bounded section
#'   "MainStreet" - whole street
#'   "Street1 [A,B], Street2, Street3 [C,D]" - multi-segment route
#' Returns list of segments, each with main and optional cross
parse_street_syntax <- function(query) {
  query <- trimws(query)
  if (nchar(query) == 0) return(list())
  
  # Split by comma OUTSIDE brackets
  # Track bracket depth to find split points
  chars <- strsplit(query, "")[[1]]
  depth <- 0
  splits <- c(0)
  
 for (i in seq_along(chars)) {
    if (chars[i] == "[") depth <- depth + 1
    if (chars[i] == "]") depth <- depth - 1
    if (chars[i] == "," && depth == 0) {
      splits <- c(splits, i)
    }
  }
  splits <- c(splits, nchar(query) + 1)
  
  # Extract segments between splits
  segments <- list()
  for (i in seq_len(length(splits) - 1)) {
    start <- splits[i] + 1
    end <- splits[i + 1] - 1
    chunk <- trimws(substr(query, start, end))
    
    if (nchar(chunk) == 0) next
    
    # Parse individual segment: "Street [Cross1, Cross2]" or just "Street"
    bracket_match <- regmatches(chunk, regexec("^(.+?)\\s*\\[(.+)\\]$", chunk))[[1]]
    
    if (length(bracket_match) == 3) {
      main <- trimws(bracket_match[2])
      cross_part <- bracket_match[3]
      cross <- trimws(strsplit(cross_part, ",")[[1]])
      segments[[length(segments) + 1]] <- list(main = main, cross = cross)
    } else {
      segments[[length(segments) + 1]] <- list(main = chunk, cross = NULL)
    }
  }
  
  segments
}

#' Find segments of main street between cross-street intersections
#' Returns dataframe with in_range flag and context segments
#' Supports: [A, B] = between A and B
#'           [A:] = from A to end
#'           [:B] = from start to B
#' Only includes segments connected to the found intersection(s)
filter_between_cross_streets <- function(net, main_street, cross_streets, tol = 5, context = 3) {
  main_q <- toupper(main_street)
  
  # Get main street segments
  main_segs <- net |>
    filter(grepl(main_q, toupper(PRI_NAME), fixed = TRUE) |
           grepl(main_q, toupper(SEC_NAME), fixed = TRUE))
  
  if (nrow(main_segs) == 0) {
    message(sprintf("No segments found for '%s'", main_street))
    main_segs$in_range <- logical(0)
    return(main_segs)
  }
  
  n_segs <- nrow(main_segs)
  
  # Handle empty or single cross street
  cross_streets <- cross_streets[nchar(trimws(cross_streets)) > 0]
  if (length(cross_streets) == 0) {
    # No bounds - order and return all
    if (n_segs > 1) {
      ord <- order_segments_by_network(main_segs$geom)
      main_segs <- main_segs[ord, ]
    }
    main_segs$in_range <- TRUE
    return(main_segs)
  }
  
  # Build adjacency for main street segments
  main_endpoints <- lapply(seq_len(n_segs), function(i) {
    coords <- wk::wk_coords(main_segs$geom[i])
    list(
      start = c(coords$x[1], coords$y[1]),
      end = c(coords$x[nrow(coords)], coords$y[nrow(coords)])
    )
  })
  
  # Build adjacency matrix
  adjacency <- vector("list", n_segs)
  for (i in seq_len(n_segs)) {
    adjacency[[i]] <- integer()
    for (j in seq_len(n_segs)) {
      if (i == j) next
      for (pi in list(main_endpoints[[i]]$start, main_endpoints[[i]]$end)) {
        for (pj in list(main_endpoints[[j]]$start, main_endpoints[[j]]$end)) {
          if (sqrt(sum((pi - pj)^2)) <= tol) {
            adjacency[[i]] <- c(adjacency[[i]], j)
            break
          }
        }
        if (j %in% adjacency[[i]]) break
      }
    }
  }
  
  # Find intersection index with a cross street
  find_intersection_idx <- function(cross_name) {
    cross_q <- toupper(trimws(cross_name))
    cross_segs <- net |>
      filter(grepl(cross_q, toupper(PRI_NAME), fixed = TRUE) |
             grepl(cross_q, toupper(SEC_NAME), fixed = TRUE))
    
    if (nrow(cross_segs) == 0) {
      message(sprintf("  Cross street '%s' not found", cross_name))
      return(NULL)
    }
    
    # Get endpoints of cross street segments
    cross_endpoints <- lapply(seq_len(nrow(cross_segs)), function(i) {
      coords <- wk::wk_coords(cross_segs$geom[i])
      list(
        start = c(coords$x[1], coords$y[1]),
        end = c(coords$x[nrow(coords)], coords$y[nrow(coords)])
      )
    })
    
    # Find shared endpoints between main and cross
    for (mi in seq_along(main_endpoints)) {
      for (mpt in list(main_endpoints[[mi]]$start, main_endpoints[[mi]]$end)) {
        for (ci in seq_along(cross_endpoints)) {
          for (cpt in list(cross_endpoints[[ci]]$start, cross_endpoints[[ci]]$end)) {
            if (sqrt(sum((mpt - cpt)^2)) <= tol) {
              return(mi)
            }
          }
        }
      }
    }
    message(sprintf("  No intersection found between '%s' and '%s'", main_street, cross_name))
    NULL
  }
  
  # Find connected component containing a given segment index
  find_connected_component <- function(start_idx) {
    visited <- logical(n_segs)
    queue <- start_idx
    component <- integer()
    
    while (length(queue) > 0) {
      current <- queue[1]
      queue <- queue[-1]
      
      if (visited[current]) next
      visited[current] <- TRUE
      component <- c(component, current)
      
      # Add unvisited neighbors
      for (neighbor in adjacency[[current]]) {
        if (!visited[neighbor]) {
          queue <- c(queue, neighbor)
        }
      }
    }
    sort(component)
  }
  
  # Find intersections
  idx1 <- find_intersection_idx(cross_streets[1])
  idx2 <- if (length(cross_streets) >= 2) find_intersection_idx(cross_streets[2]) else NULL
  
  # Get connected component from first found intersection
  anchor_idx <- if (!is.null(idx1)) idx1 else idx2
  if (is.null(anchor_idx)) {
    message(sprintf("No intersections found - showing all %s (%d segments)", main_street, n_segs))
    if (n_segs > 1) {
      ord <- order_segments_by_network(main_segs$geom)
      main_segs <- main_segs[ord, ]
    }
    main_segs$in_range <- TRUE
    return(main_segs)
  }
  
  # Filter to connected component only
  connected <- find_connected_component(anchor_idx)
  message(sprintf("  Connected component has %d of %d segments", length(connected), n_segs))
  
  main_segs <- main_segs[connected, ]
  n_segs <- nrow(main_segs)
  
  # Re-order the connected subset
  if (n_segs > 1) {
    ord <- order_segments_by_network(main_segs$geom)
    main_segs <- main_segs[ord, ]
  }
  
  # Re-build endpoints for filtered segments
  main_endpoints <- lapply(seq_len(n_segs), function(i) {
    coords <- wk::wk_coords(main_segs$geom[i])
    list(
      start = c(coords$x[1], coords$y[1]),
      end = c(coords$x[nrow(coords)], coords$y[nrow(coords)])
    )
  })
  
  # Re-find intersection indices in filtered/reordered set
  find_intersection_idx_filtered <- function(cross_name) {
    cross_q <- toupper(trimws(cross_name))
    cross_segs <- net |>
      filter(grepl(cross_q, toupper(PRI_NAME), fixed = TRUE) |
             grepl(cross_q, toupper(SEC_NAME), fixed = TRUE))
    
    if (nrow(cross_segs) == 0) return(NULL)
    
    cross_endpoints <- lapply(seq_len(nrow(cross_segs)), function(i) {
      coords <- wk::wk_coords(cross_segs$geom[i])
      list(
        start = c(coords$x[1], coords$y[1]),
        end = c(coords$x[nrow(coords)], coords$y[nrow(coords)])
      )
    })
    
    for (mi in seq_along(main_endpoints)) {
      for (mpt in list(main_endpoints[[mi]]$start, main_endpoints[[mi]]$end)) {
        for (ci in seq_along(cross_endpoints)) {
          for (cpt in list(cross_endpoints[[ci]]$start, cross_endpoints[[ci]]$end)) {
            if (sqrt(sum((mpt - cpt)^2)) <= tol) {
              return(mi)
            }
          }
        }
      }
    }
    NULL
  }
  
  # Find indices in filtered set
  idx1 <- find_intersection_idx_filtered(cross_streets[1])
  idx2 <- if (length(cross_streets) >= 2) find_intersection_idx_filtered(cross_streets[2]) else NULL
  
  # Determine range based on what we found
  if (!is.null(idx1) && !is.null(idx2)) {
    range_start <- min(idx1, idx2)
    range_end <- max(idx1, idx2)
    message(sprintf("Range: segments %d to %d", range_start, range_end))
  } else if (!is.null(idx1)) {
    if (idx1 <= n_segs / 2) {
      range_start <- 1
      range_end <- idx1
    } else {
      range_start <- idx1
      range_end <- n_segs
    }
    message(sprintf("Single bound at segment %d", idx1))
  } else if (!is.null(idx2)) {
    if (idx2 <= n_segs / 2) {
      range_start <- 1
      range_end <- idx2
    } else {
      range_start <- idx2
      range_end <- n_segs
    }
    message(sprintf("Single bound at segment %d", idx2))
  } else {
    # Shouldn't happen but fallback
    range_start <- 1
    range_end <- n_segs
  }
  
  # Add context on each side
  context_start <- max(1, range_start - context)
  context_end <- min(n_segs, range_end + context)
  
  # Flag which are in range vs context
  main_segs$in_range <- FALSE
  main_segs$in_range[range_start:range_end] <- TRUE
  
  # Return with context
  main_segs[context_start:context_end, ]
}

# --- UI ---
ui <- fluidPage(
  titlePanel("Cycle Infrastructure Mapper"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      h4("1. Search"),
      checkboxInput("hobart_only", "Hobart area only (10km)", value = TRUE),
      textInput("street_name", "Street name / route:", placeholder = "Argyle [Federal, Burnett], Weld"),
      helpText("[A, B] bounds section • [A,] from A to end • comma separates route parts"),
      checkboxInput("show_context", "Show context segments", value = TRUE),
      fluidRow(
        column(6, numericInput("anchor_lat", "Lat:", value = NA, step = 0.0001)),
        column(6, numericInput("anchor_lng", "Lng:", value = NA, step = 0.0001))
      ),
      sliderInput("buffer_m", "Search radius (m):", 100, 2000, 500, step = 100),
      actionButton("search_btn", "Search", class = "btn-primary"),
      
      hr(),
      
      h4("2. Selection"),
      verbatimTextOutput("selection_info"),
      actionButton("clear_btn", "Clear"),
      
      hr(),
      
      h4("3. Save"),
      textInput("entity_id", "Cycle path ID:"),
      textAreaInput("description", "Description:", rows = 2),
      actionButton("save_btn", "Save", class = "btn-success"),
      
      hr(),
      verbatimTextOutput("saved_summary"),
      downloadButton("export_btn", "Export CSV")
    ),
    
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "main_tabs",
        tabPanel("Map", 
          leafletOutput("map", height = "600px"),
          h4("Candidates"),
          DT::dataTableOutput("candidates_tbl")
        ),
        tabPanel("Progress",
          leafletOutput("progress_map", height = "600px"),
          h4("Saved Mappings"),
          DT::dataTableOutput("progress_tbl")
        )
      )
    )
  )
)

# --- Server ---
server <- function(input, output, session) {
  
  # Load network once
  network <- reactiveVal(NULL)
  
  observe({
    showModal(modalDialog("Loading network...", footer = NULL))
    net <- tryCatch(load_network(), error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
      NULL
    })
    network(net)
    removeModal()
  })
  
  # Search results
  candidates <- reactiveVal(data.frame())
  
  # Selected TRANSEG_IDs
  selected <- reactiveVal(character())
  
  # Saved mappings - load from file on startup
  mappings <- reactiveVal(load_mappings())
  
  # Search
  observeEvent(input$search_btn, {
    req(network())
    net <- network()
    
    results <- net
    
    # Default Hobart area filter (10km radius)
    if (isTRUE(input$hobart_only)) {
      hobart_buf <- 10000 / 111000  # ~10km in degrees
      results <- results |>
        filter(abs(clng - HOBART_CENTER[2]) < hobart_buf,
               abs(clat - HOBART_CENTER[1]) < hobart_buf)
    }
    
    # Parse search syntax - now returns list of segments
    parsed_segments <- parse_street_syntax(input$street_name)
    
    if (length(parsed_segments) > 0) {
      # Process each segment and combine
      all_results <- list()
      
      for (seg in parsed_segments) {
        if (!is.null(seg$cross) && length(seg$cross) >= 2) {
          # Bounded segment
          seg_result <- filter_between_cross_streets(results, seg$main, seg$cross, context = 3)
        } else {
          # Whole street
          q <- toupper(seg$main)
          seg_result <- results |>
            filter(grepl(q, toupper(PRI_NAME), fixed = TRUE) |
                   grepl(q, toupper(SEC_NAME), fixed = TRUE))
          
          # Order it
          if (nrow(seg_result) > 1) {
            ord <- order_segments_by_network(seg_result$geom)
            seg_result <- seg_result[ord, ]
          }
          seg_result$in_range <- TRUE
        }
        
        if (nrow(seg_result) > 0) {
          all_results[[length(all_results) + 1]] <- seg_result
        }
      }
      
      # Combine all segments - use bind_rows for arrow compatibility
      if (length(all_results) > 0) {
        results <- dplyr::bind_rows(all_results)
        # Remove duplicates (keep first occurrence)
        results <- results[!duplicated(results$TRANSEG_ID), ]
      } else {
        results <- results[0, ]  # Empty
      }
      
      showNotification(sprintf("Found %d segments across %d route parts", 
                               nrow(results), length(parsed_segments)), 
                       duration = 3)
    } else {
      results$in_range <- TRUE
    }
    
    # Filter by proximity to anchor point
    if (!is.na(input$anchor_lng) && !is.na(input$anchor_lat)) {
      buf <- input$buffer_m / 111000  # rough degrees
      results <- results |>
        filter(abs(clng - input$anchor_lng) < buf,
               abs(clat - input$anchor_lat) < buf)
    }
    
    showNotification(sprintf("Found %d segments", nrow(results)), duration = 2)
    
    # Order by network connectivity - but only for single segments
    # Multi-segment routes preserve the order from the query
    if (nrow(results) > 1 && length(parsed_segments) <= 1) {
      ord <- order_segments_by_network(results$geom)
      results <- results[ord, ]
    }
    results$seq_order <- seq_len(nrow(results))
    
    candidates(results)
    selected(character())
    
    # Zoom to results
    if (nrow(results) > 0) {
      leafletProxy("map") |>
        fitBounds(min(results$clng) - 0.002, min(results$clat) - 0.002,
                  max(results$clng) + 0.002, max(results$clat) + 0.002)
    }
  })
  
  # Clear
  observeEvent(input$clear_btn, {
    selected(character())
  })
  
  # Save - updates in-memory and writes to file
  observeEvent(input$save_btn, {
    req(input$entity_id, length(selected()) > 0)
    
    m <- mappings()
    m[[input$entity_id]] <- list(
      id = input$entity_id,
      description = input$description,
      transeg_ids = selected()
    )
    mappings(m)
    
    # Auto-save to file
    save_mappings(m)
    
    showNotification(sprintf("Saved %s with %d segments", 
                             input$entity_id, length(selected())))
  })
  
  # Selection info with total length
  output$selection_info <- renderText({
    sel <- selected()
    cands <- candidates()
    if (length(sel) == 0) return("None selected")
    
    # Calculate total length
    total_len <- sum(cands$COMP_LEN[as.character(cands$TRANSEG_ID) %in% sel], na.rm = TRUE)
    
    sprintf("%d segments (%.0fm total):\n%s%s", 
            length(sel), 
            total_len,
            paste(head(sel, 5), collapse = "\n"),
            if (length(sel) > 5) "\n..." else "")
  })
  
  # Saved summary
  output$saved_summary <- renderText({
    m <- mappings()
    if (length(m) == 0) return("(none)")
    paste(sapply(m, function(x) sprintf("%s: %d", x$id, length(x$transeg_ids))), 
          collapse = "\n")
  })
  
  # Export
  output$export_btn <- downloadHandler(
    filename = function() paste0("cycle_mappings_", Sys.Date(), ".csv"),
    content = function(file) {
      m <- mappings()
      if (length(m) == 0) {
        write.csv(data.frame(id = character(), transeg_id = character()), file)
        return()
      }
      out <- do.call(rbind, lapply(m, function(x) {
        data.frame(id = x$id, description = x$description, 
                   transeg_id = x$transeg_ids, stringsAsFactors = FALSE)
      }))
      write.csv(out, file, row.names = FALSE)
    }
  )
  
  # Base map
  output$map <- renderLeaflet({
    leaflet() |>
      addProviderTiles(providers$CartoDB.Positron) |>
      setView(HOBART_CENTER[2], HOBART_CENTER[1], zoom = 14)
  })
  
  # Draw candidates and selected - reactive to both and show_context
  observe({
    cands <- candidates()
    sel <- selected()
    show_ctx <- isTRUE(input$show_context)
    
    # Filter to in_range only if show_context is off
    if (!show_ctx && "in_range" %in% names(cands) && nrow(cands) > 0) {
      # Keep only in_range, but don't lose selected that might be in context
      keep <- cands$in_range | (as.character(cands$TRANSEG_ID) %in% sel)
      cands <- cands[keep, ]
    }
    
    # Debug
    message(sprintf("Drawing: %d candidates, %d selected, show_context=%s", 
                    nrow(cands), length(sel), show_ctx))
    
    proxy <- leafletProxy("map") |>
      clearGroup("candidates") |>
      clearGroup("context") |>
      clearGroup("selected") |>
      clearGroup("anchor")
    
    if (nrow(cands) > 0) {
      has_in_range <- "in_range" %in% names(cands)
      
      # Draw context segments first (lightest, underneath)
      if (has_in_range && show_ctx) {
        ctx_idx <- which(!cands$in_range)
        for (i in ctx_idx) {
          id_str <- as.character(cands$TRANSEG_ID[i])
          ll <- geom_to_lonlat(cands$geom[i])
          seq_num <- if ("seq_order" %in% names(cands)) cands$seq_order[i] else i
          
          proxy <- proxy |>
            addPolylines(
              lng = ll$lng, lat = ll$lat,
              color = STYLE$context$color,
              weight = STYLE$context$weight,
              opacity = STYLE$context$opacity,
              group = "context",
              layerId = paste0("cand_", id_str),
              popup = sprintf("<b>#%d %s</b> (context)<br>ID: %s<br>Length: %.0fm", 
                             seq_num, cands$PRI_NAME[i], id_str, cands$COMP_LEN[i])
            )
        }
      }
      
      # Draw in-range candidates (grey, middle layer)
      in_range_idx <- if (has_in_range) which(cands$in_range) else seq_len(nrow(cands))
      for (i in in_range_idx) {
        id_str <- as.character(cands$TRANSEG_ID[i])
        if (id_str %in% sel) next  # Skip selected, draw on top
        
        ll <- geom_to_lonlat(cands$geom[i])
        seq_num <- if ("seq_order" %in% names(cands)) cands$seq_order[i] else i
        
        proxy <- proxy |>
          addPolylines(
            lng = ll$lng, lat = ll$lat,
            color = STYLE$candidate$color,
            weight = STYLE$candidate$weight,
            opacity = STYLE$candidate$opacity,
            group = "candidates",
            layerId = paste0("cand_", id_str),
            popup = sprintf("<b>#%d %s</b><br>ID: %s<br>Length: %.0fm", 
                           seq_num, cands$PRI_NAME[i], id_str, cands$COMP_LEN[i])
          )
      }
      
      # Draw selected ON TOP (bright green, thicker)
      sel_idx <- which(as.character(cands$TRANSEG_ID) %in% sel)
      if (length(sel_idx) > 0) {
        for (i in sel_idx) {
          id_str <- as.character(cands$TRANSEG_ID[i])
          ll <- geom_to_lonlat(cands$geom[i])
          seq_num <- if ("seq_order" %in% names(cands)) cands$seq_order[i] else i
          
          proxy <- proxy |>
            addPolylines(
              lng = ll$lng, lat = ll$lat,
              color = STYLE$selected$color,
              weight = STYLE$selected$weight,
              opacity = STYLE$selected$opacity,
              group = "selected",
              layerId = paste0("sel_", id_str),
              popup = sprintf("<b>#%d %s</b><br>ID: %s<br>Length: %.0fm<br><b>✓ SELECTED</b>", 
                             seq_num, cands$PRI_NAME[i], id_str, cands$COMP_LEN[i])
            )
        }
      }
    }
    
    # Anchor point
    if (!is.na(input$anchor_lng) && !is.na(input$anchor_lat)) {
      proxy |> addCircleMarkers(
        input$anchor_lng, input$anchor_lat,
        radius = STYLE$anchor$radius, 
        color = STYLE$anchor$color, 
        fillOpacity = 0.8,
        group = "anchor",
        popup = "Search anchor"
      )
    }
  })
  
  # Click to select/deselect with feedback
  observeEvent(input$map_shape_click, {
    click_id <- input$map_shape_click$id
    if (!is.null(click_id)) {
      # Strip prefix (cand_ or sel_) to get actual TRANSEG_ID
      id <- sub("^(cand_|sel_)", "", click_id)
      
      sel <- selected()
      if (id %in% sel) {
        selected(setdiff(sel, id))
        showNotification(sprintf("Deselected %s", id), type = "warning", duration = 2)
      } else {
        selected(c(sel, id))
        showNotification(sprintf("Selected %s", id), type = "message", duration = 2)
      }
    }
  })
  
  # Candidates table - scrollable, auto-sized, respects context toggle
  # NOTE: Don't depend on selected() here or it creates circular update
  output$candidates_tbl <- DT::renderDataTable({
    cands <- candidates()
    show_ctx <- isTRUE(input$show_context)
    
    if (nrow(cands) == 0) return(NULL)
    
    has_in_range <- "in_range" %in% names(cands)
    
    # Filter if not showing context (but we can't filter by selected here)
    if (!show_ctx && has_in_range) {
      cands <- cands[cands$in_range, ]
    }
    
    # Build display table
    display_cols <- c("seq_order", "TRANSEG_ID", "PRI_NAME", "SEC_NAME", "COMP_LEN")
    if (has_in_range && show_ctx) {
      display_cols <- c(display_cols, "in_range")
    }
    
    tbl <- cands |> select(any_of(display_cols))
    
    dt <- DT::datatable(
      tbl,
      selection = "multiple", 
      options = list(
        pageLength = max(nrow(tbl), 10),
        scrollX = TRUE,
        scrollY = "300px",
        order = list(list(0, 'asc'))
      )
    )
    
    # Only format in_range if it exists
    if (has_in_range && "in_range" %in% names(tbl)) {
      dt <- dt |>
        DT::formatStyle(
          "in_range",
          backgroundColor = DT::styleEqual(c(TRUE, FALSE), c("#e8f5e9", "#f5f5f5"))
        )
    }
    
    dt
  })
  
  # Table selection syncs to map - must match filtered table
  observeEvent(input$candidates_tbl_rows_selected, {
    cands <- candidates()
    show_ctx <- isTRUE(input$show_context)
    rows <- input$candidates_tbl_rows_selected
    
    if (nrow(cands) == 0) return()
    
    has_in_range <- "in_range" %in% names(cands)
    
    # Apply same filter as table (without selected dependency)
    if (!show_ctx && has_in_range) {
      cands <- cands[cands$in_range, ]
    }
    
    if (length(rows) > 0 && max(rows) <= nrow(cands)) {
      selected(as.character(cands$TRANSEG_ID[rows]))
    } else if (length(rows) == 0) {
      # Don't clear selection when rows is empty - that happens on re-render
    }
  }, ignoreNULL = FALSE)
  
  # --- Progress tab ---
  
  # Progress map
  output$progress_map <- renderLeaflet({
    leaflet() |>
      addProviderTiles(providers$CartoDB.Positron) |>
      setView(HOBART_CENTER[2], HOBART_CENTER[1], zoom = 13)
  })
  
  # Draw saved mappings on progress map - trigger on tab switch too
  observe({
    # Trigger when switching to Progress tab
    req(input$main_tabs == "Progress")
    
    m <- mappings()
    net <- network()
    req(net, length(m) > 0)
    
    message(sprintf("Progress map: drawing %d mappings", length(m)))
    
    proxy <- leafletProxy("progress_map") |>
      clearGroup("saved")
    
    # Get all saved transeg_ids with their mapping id
    all_saved <- do.call(rbind, lapply(m, function(x) {
      data.frame(map_id = x$id, transeg_id = as.character(x$transeg_ids), stringsAsFactors = FALSE)
    }))
    
    message(sprintf("  Total segments to draw: %d", nrow(all_saved)))
    
    # Color palette for different mapping IDs
    unique_ids <- unique(all_saved$map_id)
    pal <- colorFactor(palette = "Set1", domain = unique_ids)
    
    # Collect bounds for zoom
    all_lngs <- c()
    all_lats <- c()
    
    # Draw each saved segment
    drawn <- 0
    for (i in seq_len(nrow(all_saved))) {
      tid <- all_saved$transeg_id[i]
      mid <- all_saved$map_id[i]
      
      seg_row <- which(as.character(net$TRANSEG_ID) == tid)
      if (length(seg_row) == 0) {
        message(sprintf("  TRANSEG %s not found in network!", tid))
        next
      }
      
      ll <- geom_to_lonlat(net$geom[seg_row[1]])
      all_lngs <- c(all_lngs, ll$lng)
      all_lats <- c(all_lats, ll$lat)
      
      proxy <- proxy |>
        addPolylines(
          lng = ll$lng, lat = ll$lat,
          color = pal(mid),
          weight = STYLE$saved$weight,
          opacity = STYLE$saved$opacity,
          group = "saved",
          popup = sprintf("<b>%s</b><br>TRANSEG: %s", mid, tid)
        )
      drawn <- drawn + 1
    }
    message(sprintf("  Drew %d segments", drawn))
    
    # Zoom to fit all saved segments
    if (length(all_lngs) > 0) {
      proxy |> fitBounds(min(all_lngs) - 0.002, min(all_lats) - 0.002,
                         max(all_lngs) + 0.002, max(all_lats) + 0.002)
    }
  })
  
  # Progress table
  output$progress_tbl <- DT::renderDataTable({
    m <- mappings()
    if (length(m) == 0) return(NULL)
    
    # Summary by mapping ID
    summary_df <- do.call(rbind, lapply(m, function(x) {
      data.frame(
        id = x$id,
        description = x$description,
        n_segments = length(x$transeg_ids),
        stringsAsFactors = FALSE
      )
    }))
    
    DT::datatable(summary_df, options = list(pageLength = 10))
  })
}


shinyApp(ui, server)
