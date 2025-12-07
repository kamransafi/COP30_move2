#' Download and annotate Copernicus DSM (30 m) for a move2 object
#'
#' This is the main user-facing function. It automatically decides whether to
#' request the DSM directly (for small extents) or use tiling (for large extents)
#' based on the requested area size and a pixel-count threshold.
#'
#' For small spatial extents this function delegates directly to the single-request
#' implementation. For large extents it internally splits the global bounding box
#' into regular tiles, requests a DSM per tile, and merges the per-location
#' elevation values back into a single move2 object.
#'
#' Three tiling strategies are available:
#' \itemize{
#'   \item \code{"grid"} (default): Regular grid over entire bounding box.
#'   \item \code{"data-aware"}: Only tiles containing movement data are requested.
#'   \item \code{"convex-hull"}: Only tiles intersecting the convex hull of points.
#' }
#'
#' Authentication is done via the \pkg{CDSE} package using a Sentinel Hub
#' OAuth client ID and secret stored securely in the system keyring via
#' the \pkg{keyring} package (entries named \code{"sentinelhub_client_id"}
#' and \code{"sentinelhub_client_secret"} under the chosen \code{key_service}).
#'
#' @param x A \code{move2} object containing animal trajectory data.
#' @param extent_factor Numeric > 0. Factor by which the bounding box of \code{x}
#'   should be expanded before requesting the DSM (1 = original extent, 1.1 = +10\%).
#' @param method Character; interpolation method used by \code{terra::extract()}.
#'   One of \code{"nearest"} (nearest neighbour; internally \code{"simple"}) or
#'   \code{"bilinear"} (bilinear interpolation). Default: \code{"bilinear"}.
#' @param max_pixels_per_request Numeric; approximate maximum number of
#'   DSM pixels per Process API request before tiling is triggered.
#'   Default: 5e6 pixels (~0.5° × 0.5° at 30 m resolution).
#' @param tile_dx,tile_dy Numeric; tile width/height in degrees (lon/lat)
#'   used when tiling is required. Defaults: 0.25 degrees each (~28 km at equator).
#' @param tiling_strategy Character; which tiles to request when tiling is used.
#'   One of \code{"grid"} (all tiles in bbox), \code{"data-aware"} (only tiles
#'   with points), or \code{"convex-hull"} (tiles intersecting convex hull).
#'   Default: \code{"grid"}.
#' @param buffer_deg Numeric; buffer in degrees around points/tiles when using
#'   "data-aware" or "convex-hull" strategies to ensure interpolation doesn't
#'   fail at tile boundaries. Default: 0.01 degrees (~1 km at equator).
#'   Increase this if you experience NA values near tile edges.
#' @param plot_overview Logical; if \code{TRUE} (default), produce a base R
#'   overview plot with the DSM as background and \code{move2} points coloured
#'   by track ID. For tiled requests, one final plot is drawn on the merged result.
#' @param dem_file Optional character file path for the output GeoTIFF.
#'   If \code{NULL} (default), a temporary file is created. Only used in single-request mode.
#' @param key_service Character; name of the keyring service under which the
#'   Sentinel Hub client ID and secret are stored. Required; no default.
#'
#' @return An (invisible) list with:
#'   \itemize{
#'     \item \code{move2}: the input \code{move2} object with a new numeric column
#'           \code{dem_copernicus30} containing DSM elevations (in metres) at
#'           each location (NA where geometry is missing or outside DEM).
#'     \item \code{dem}: a \code{terra::SpatRaster} containing the downloaded DSM
#'           (for single-request mode) or the first tile (for tiled mode).
#'     \item \code{dem_list}: list of \code{terra::SpatRaster} tiles (NULL for single-request).
#'     \item \code{demfile}: character path to the GeoTIFF file on disk (single-request mode).
#'     \item \code{demfiles}: character vector of tile GeoTIFF paths (tiled mode).
#'     \item \code{tiles}: data frame of tile extents (x0, y0, x1, y1); single row for single-request.
#'     \item \code{tiled}: logical; TRUE if tiling was used, FALSE otherwise.
#'     \item \code{tiling_strategy}: character; strategy used (if tiled).
#'   }
#'
#' @details
#' The function assumes that the DEM is requested and returned in WGS84
#' longitude/latitude (EPSG:4326). If you change the DEM CRS in the Process API request,
#' you must adapt the extraction step accordingly.
#'
#' This function is intended as a convenience wrapper for exploratory analysis
#' and not as a production-grade client. Use with appropriate API quotas
#' and follow Copernicus Data Space terms of use.
#'
#' @examples
#' \dontrun{
#' library(move2)
#'
#' # Example with bundled move2 data (small extent, single request)
#' fishers <- mt_read(mt_example())
#'
#' res <- get_COP30_move2(
#'   x             = fishers,
#'   extent_factor = 1.2,
#'   method        = "bilinear",
#'   plot_overview = TRUE,
#'   key_service   = "SentinelHUB"
#' )
#'
#' fishers_annotated <- res$move2
#' head(fishers_annotated$dem_copernicus30)
#'
#' # For larger extents with long-distance migration, use data-aware tiling:
#' # res <- get_COP30_move2(large_migratory_object, 
#' #                        tiling_strategy = "data-aware",
#' #                        key_service = "SentinelHUB")
#' }
#'
#' @export
get_COP30_move2 <- function(x,
                            extent_factor          = 1.1,
                            method                 = c("nearest", "bilinear"),
                            plot_overview          = TRUE,
                            max_pixels_per_request = 5e6,
                            tile_dx                = 0.25,
                            tile_dy                = 0.25,
                            tiling_strategy        = c("grid", "data-aware", "convex-hull"),
                            buffer_deg             = 0.01,
                            dem_file               = NULL,
                            key_service            = NULL) {
  
  method <- match.arg(method)
  tiling_strategy <- match.arg(tiling_strategy)
  
  if (!inherits(x, "move2")) {
    stop("x must be a move2 object.")
  }
  
  if (is.null(key_service) || !is.character(key_service) || key_service == "") {
    stop("key_service must be a non-empty character string naming a keyring service.")
  }
  
  if (buffer_deg < 0) {
    stop("buffer_deg must be non-negative.")
  }
  
  # Compute global bounding box in WGS84
  x_wgs <- sf::st_transform(x, 4326)
  global_bbox <- .compute_global_bbox(x_wgs, extent_factor)
  
  # Estimate pixel count and decide whether tiling is needed
  n_pix <- .estimate_dsm_pixels(global_bbox)
  
  if (is.na(n_pix) || n_pix <= max_pixels_per_request) {
    # --- Single-request path ---
    res <- .get_dsm_single_request(
      x             = x,
      bbox_override = global_bbox,
      method        = method,
      plot_overview = plot_overview,
      dem_file      = dem_file,
      key_service   = key_service
    )
    
    return(invisible(list(
      move2            = res$move2,
      dem              = res$dem,
      dem_list         = NULL,
      demfile          = res$demfile,
      demfiles         = NULL,
      tiles            = data.frame(
        x0 = global_bbox[1],
        y0 = global_bbox[2],
        x1 = global_bbox[3],
        y1 = global_bbox[4]
      ),
      tiled            = FALSE,
      tiling_strategy  = NA_character_
    )))
  }
  
  # --- Tiling path ---
  res <- .get_dsm_tiled(
    x               = x,
    global_bbox     = global_bbox,
    method          = method,
    plot_overview   = plot_overview,
    tile_dx         = tile_dx,
    tile_dy         = tile_dy,
    tiling_strategy = tiling_strategy,
    buffer_deg      = buffer_deg,
    key_service     = key_service
  )
  
  invisible(list(
    move2            = res$move2,
    dem              = res$dem_list[[1]],
    dem_list         = res$dem_list,
    demfile          = NULL,
    demfiles         = res$demfiles,
    tiles            = res$tiles,
    tiled            = TRUE,
    tiling_strategy  = tiling_strategy
  ))
}


# =========================================================================
# Internal: Tiling strategy helpers
# =========================================================================

#' Internal: Filter tiles by strategy with buffer zones
#'
#' @param tiles Data frame of all potential tiles.
#' @param sf_pts sf object with point locations.
#' @param strategy Character; "grid", "data-aware", or "convex-hull".
#' @param buffer_deg Numeric; buffer in degrees around points to ensure
#'   interpolation doesn't fail at tile edges. Default 0.01 degrees (~1 km).
#'
#' @return Data frame of tiles to request.
#'
#' @keywords internal
.filter_tiles_by_strategy <- function(tiles, sf_pts, strategy = "grid", 
                                      buffer_deg = 0.01) {
  
  if (strategy == "grid") {
    # Request all tiles
    return(tiles)
  }
  
  if (strategy == "data-aware") {
    # Only tiles that have points or their buffer zones
    tiles_with_data <- integer()
    
    for (i in seq_len(nrow(tiles))) {
      tile_bbox <- as.numeric(tiles[i, c("x0", "y0", "x1", "y1")])
      
      # Expand tile by buffer to catch points near the edge
      poly_buffered <- sf::st_as_sfc(sf::st_bbox(c(
        xmin = tile_bbox[1] - buffer_deg,
        ymin = tile_bbox[2] - buffer_deg,
        xmax = tile_bbox[3] + buffer_deg,
        ymax = tile_bbox[4] + buffer_deg
      ), crs = 4326))
      
      if (any(lengths(sf::st_within(sf_pts, poly_buffered)) > 0)) {
        tiles_with_data <- c(tiles_with_data, i)
      }
    }
    
    return(tiles[tiles_with_data, ])
  }
  
  if (strategy == "convex-hull") {
    # Only tiles intersecting the buffered convex hull
    hull <- sf::st_convex_hull(sf::st_union(sf::st_geometry(sf_pts)))
    hull_buffered <- sf::st_buffer(hull, dist = buffer_deg)
    tiles_with_hull <- integer()
    
    for (i in seq_len(nrow(tiles))) {
      tile_bbox <- as.numeric(tiles[i, c("x0", "y0", "x1", "y1")])
      poly <- sf::st_as_sfc(sf::st_bbox(c(
        xmin = tile_bbox[1],
        ymin = tile_bbox[2],
        xmax = tile_bbox[3],
        ymax = tile_bbox[4]
      ), crs = 4326))
      
      if (!sf::st_is_empty(sf::st_intersection(poly, hull_buffered))) {
        tiles_with_hull <- c(tiles_with_hull, i)
      }
    }
    
    return(tiles[tiles_with_hull, ])
  }
  
  stop("Unknown tiling_strategy: ", strategy)
}


# =========================================================================
# Internal: Common operations
# =========================================================================

#' Internal: Compute global bounding box
#'
#' @param x_wgs move2 object in EPSG:4326.
#' @param extent_factor Expansion factor.
#'
#' @return Numeric vector c(xmin, ymin, xmax, ymax).
#'
#' @keywords internal
.compute_global_bbox <- function(x_wgs, extent_factor) {
  bb <- sf::st_bbox(x_wgs)
  cx <- (bb["xmin"] + bb["xmax"]) / 2
  cy <- (bb["ymin"] + bb["ymax"]) / 2
  dx <- (bb["xmax"] - bb["xmin"]) * extent_factor / 2
  dy <- (bb["ymax"] - bb["ymin"]) * extent_factor / 2
  c(cx - dx, cy - dy, cx + dx, cy + dy)
}


#' Internal: Retrieve Sentinel Hub OAuth token
#'
#' @param key_service Character; keyring service name.
#'
#' @return Character; OAuth token.
#'
#' @keywords internal
.get_sentinelhub_token <- function(key_service) {
  sh_client_id  <- keyring::key_get("sentinelhub_client_id",
                                    service = key_service)
  sh_client_sec <- keyring::key_get("sentinelhub_client_secret",
                                    service = key_service)
  
  if (is.null(sh_client_id) || is.null(sh_client_sec)) {
    stop("Sentinel Hub credentials not found in keyring under service '",
         key_service, "'.")
  }
  
  CDSE::GetOAuthToken(id = sh_client_id, secret = sh_client_sec)
}


#' Internal: Request DSM from Process API
#'
#' Sends a Process API request to Sentinel Hub for a given bounding box
#' and returns the raw TIFF response.
#'
#' @param bbox_vec Numeric vector c(xmin, ymin, xmax, ymax).
#' @param token Character; OAuth token.
#'
#' @return Raw response body (TIFF bytes).
#'
#' @keywords internal
.call_process_api <- function(bbox_vec, token) {
  process_body <- list(
    input = list(
      bounds = list(
        properties = list(
          crs = "http://www.opengis.net/def/crs/OGC/1.3/CRS84"
        ),
        bbox = as.numeric(bbox_vec)
      ),
      data = list(list(
        type = "dem",
        dataFilter = list(
          demInstance = "COPERNICUS_30"
        ),
        processing = list(
          upsampling   = "BILINEAR",
          downsampling = "BILINEAR"
        )
      ))
    ),
    output = list(
      resx = 0.0003,
      resy = 0.0003,
      responses = list(list(
        identifier = "default",
        format = list(
          type = "image/tiff"
        )
      ))
    ),
    evalscript = .dem_evalscript
  )
  
  process_url <- "https://sh.dataspace.copernicus.eu/api/v1/process"
  
  resp <- httr2::request(process_url) |>
    httr2::req_headers(
      Authorization  = paste("Bearer", token),
      `Content-Type` = "application/json"
    ) |>
    httr2::req_body_json(process_body, auto_unbox = TRUE) |>
    httr2::req_perform()
  
  httr2::resp_check_status(resp)
  httr2::resp_body_raw(resp)
}


#' Internal: Extract DEM values at locations
#'
#' @param x_wgs move2 object in EPSG:4326.
#' @param dem_rast SpatRaster with DEM.
#' @param method Character; "nearest" or "bilinear".
#'
#' @return Numeric vector of DEM values (with NAs for invalid points).
#'
#' @keywords internal
.extract_dem_values <- function(x_wgs, dem_rast, method) {
  sf_pts <- sf::st_as_sf(x_wgs)
  coords <- sf::st_coordinates(sf_pts)
  valid  <- !is.na(coords[, "X"]) & !is.na(coords[, "Y"])
  
  dem_vals <- rep(NA_real_, nrow(sf_pts))
  
  if (any(valid)) {
    pts_valid <- terra::vect(sf_pts[valid, ])
    method_terra <- if (method == "nearest") "simple" else "bilinear"
    dem_valid    <- terra::extract(dem_rast, pts_valid,
                                   method = method_terra)[, 2]
    dem_vals[valid] <- dem_valid
  }
  
  dem_vals
}


# =========================================================================
# Internal: Single-request DSM retrieval
# =========================================================================

#' Internal: Download DSM for a single bounding box
#'
#' @param x A \code{move2} object.
#' @param bbox_override Numeric vector c(xmin, ymin, xmax, ymax) in EPSG:4326.
#' @param method Character; "nearest" or "bilinear".
#' @param plot_overview Logical; if TRUE, plot the result.
#' @param dem_file Optional character file path for the output GeoTIFF.
#' @param key_service Character; keyring service name.
#'
#' @return List with move2, dem, demfile.
#'
#' @keywords internal
.get_dsm_single_request <- function(x,
                                    bbox_override,
                                    method,
                                    plot_overview,
                                    dem_file,
                                    key_service) {
  
  x_wgs <- sf::st_transform(x, 4326)
  
  # Get token and call API
  token <- .get_sentinelhub_token(key_service)
  dem_tiff <- .call_process_api(bbox_override, token)
  
  # Save and load raster
  if (is.null(dem_file)) {
    dem_file <- tempfile(fileext = ".tif")
  }
  writeBin(dem_tiff, dem_file)
  dem_rast <- terra::rast(dem_file)
  
  # Extract values and annotate
  dem_vals <- .extract_dem_values(x_wgs, dem_rast, method)
  x$dem_copernicus30 <- dem_vals
  
  # Optional plot
  if (plot_overview) {
    .plot_move2_dsm(x_wgs = x_wgs, dem_rast = dem_rast)
  }
  
  list(
    move2   = x,
    dem     = dem_rast,
    demfile = dem_file
  )
}


# =========================================================================
# Internal: Tiled DSM retrieval
# =========================================================================

#' Internal: Download DSM via tiling
#'
#' @param x A \code{move2} object.
#' @param global_bbox Numeric vector c(xmin, ymin, xmax, ymax) in EPSG:4326.
#' @param method Character; "nearest" or "bilinear".
#' @param plot_overview Logical; if TRUE, plot the merged result.
#' @param tile_dx,tile_dy Numeric; tile size in degrees.
#' @param tiling_strategy Character; strategy to use.
#' @param buffer_deg Numeric; buffer in degrees.
#' @param key_service Character; keyring service name.
#'
#' @return List with move2, dem_list, demfiles, tiles.
#'
#' @keywords internal
.get_dsm_tiled <- function(x,
                           global_bbox,
                           method,
                           plot_overview,
                           tile_dx,
                           tile_dy,
                           tiling_strategy,
                           buffer_deg,
                           key_service) {
  
  x_wgs <- sf::st_transform(x, 4326)
  all_tiles <- .make_tiles(global_bbox, tile_dx = tile_dx, tile_dy = tile_dy)
  sf_pts <- sf::st_as_sf(x_wgs)
  
  # Filter tiles based on strategy with buffer
  tiles <- .filter_tiles_by_strategy(all_tiles, sf_pts, 
                                     strategy = tiling_strategy,
                                     buffer_deg = buffer_deg)
  
  x_out <- x
  x_out$dem_copernicus30 <- NA_real_
  
  dem_list  <- vector("list", nrow(tiles))
  demfiles  <- character(nrow(tiles))
  
  message(sprintf("Tiling strategy '%s' (buffer=%.4f°): requesting %d of %d tiles",
                  tiling_strategy, buffer_deg, nrow(tiles), nrow(all_tiles)))
  
  # Process each selected tile
  for (i in seq_len(nrow(tiles))) {
    tile_bbox <- as.numeric(tiles[i, c("x0", "y0", "x1", "y1")])
    
    # Find points within this tile (with buffer for safety)
    poly <- sf::st_as_sfc(sf::st_bbox(c(
      xmin = tile_bbox[1] - buffer_deg,
      ymin = tile_bbox[2] - buffer_deg,
      xmax = tile_bbox[3] + buffer_deg,
      ymax = tile_bbox[4] + buffer_deg
    ), crs = 4326))
    
    in_tile <- lengths(sf::st_within(sf_pts, poly)) > 0
    
    if (!any(in_tile)) next
    
    x_sub <- x[in_tile, ]
    
    # Request DSM for this tile (without buffer—API gets exact tile)
    tile_bbox_exact <- as.numeric(tiles[i, c("x0", "y0", "x1", "y1")])
    res_tile <- .get_dsm_single_request(
      x             = x_sub,
      bbox_override = tile_bbox_exact,
      method        = method,
      plot_overview = FALSE,
      dem_file      = NULL,
      key_service   = key_service
    )
    
    # Merge results
    x_out$dem_copernicus30[in_tile] <- res_tile$move2$dem_copernicus30
    dem_list[[i]] <- res_tile$dem
    demfiles[i]   <- res_tile$demfile
  }
  
  # Optional final plot
  if (plot_overview) {
    i_first <- which(!vapply(dem_list, is.null, logical(1)))[1]
    if (length(i_first) > 0L) {
      dem_plot <- dem_list[[i_first]]
      .plot_move2_dsm(x_wgs = x_wgs, dem_rast = dem_plot,
                      main = "Copernicus DSM (GLO-30) with tiled tracks")
    }
  }
  
  list(
    move2    = x_out,
    dem_list = dem_list,
    demfiles = demfiles,
    tiles    = tiles
  )
}


# =========================================================================
# Internal: Plotting helper
# =========================================================================

#' Internal: Quick overview plot with split panels
#'
#' @param x_wgs A \code{move2} object in EPSG:4326.
#' @param dem_rast A \code{terra::SpatRaster} with DSM.
#' @param main Character; plot title.
#'
#' @keywords internal
.plot_move2_dsm <- function(x_wgs, dem_rast,
                            main = "Copernicus DSM (GLO-30) with tracks") {
  
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  
  # Create layout: main plot takes 80% width, legend panel takes 20%
  layout(matrix(c(1, 2), nrow = 1, ncol = 2),
         widths = c(0.8, 0.2),
         heights = c(1))
  
  # Panel 1: DSM + tracks
  plot(dem_rast, main = main, col = terrain.colors(64))
  
  track_ids <- as.factor(move2::mt_track_id(x_wgs))
  track_colors <- rainbow(length(levels(track_ids)))
  cols <- track_colors[track_ids]
  
  plot(sf::st_geometry(x_wgs), add = TRUE, col = cols, pch = 16, cex = 0.6)
  
  # Panel 2: Manual legend with ID then colored line (centered vertically)
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  plot(0, 0, type = "n", xlim = 0:1, ylim = 0:1, 
       axes = FALSE, xlab = "", ylab = "")
  
  # Sort track IDs alphabetically
  all_ids <- levels(track_ids)
  sorted_ids <- sort(all_ids)
  n_tracks <- length(sorted_ids)
  line_spacing <- 0.04
  
  # Center legend vertically: calculate total height needed and center it
  total_height <- (n_tracks - 1) * line_spacing
  center_y <- 0.5
  start_y <- center_y + total_height / 2
  
  y_pos <- start_y - seq(0, n_tracks - 1) * line_spacing
  
  for (j in seq_along(sorted_ids)) {
    track_name <- sorted_ids[j]
    # Find the index of this track in the original levels
    track_idx <- which(all_ids == track_name)
    
    # Draw track ID text (left side, right-aligned)
    text(0.45, y_pos[j], track_name, 
         cex = 0.9, adj = 1, font = 1, col = "black")
    
    # Draw colored line (right side)
    lines(c(0.50, 0.75), c(y_pos[j], y_pos[j]), 
          col = track_colors[track_idx], lwd = 4)
  }
  
  invisible(NULL)
}


# =========================================================================
# Helper functions for tiling
# =========================================================================

#' Internal: Estimate pixel count for a DSM request
#'
#' @param bbox_vec Numeric vector c(xmin, ymin, xmax, ymax) in degrees.
#' @param resx,resy Numeric; resolution in degrees (default 0.0003 ≈ 30 m).
#'
#' @return Numeric; approximate total pixel count.
#'
#' @keywords internal
.estimate_dsm_pixels <- function(bbox_vec, resx = 0.0003, resy = 0.0003) {
  nx <- (bbox_vec[3] - bbox_vec[1]) / resx
  ny <- (bbox_vec[4] - bbox_vec[2]) / resy
  as.numeric(nx * ny)
}


#' Internal: Create a regular grid of tiles
#'
#' @param bbox_vec Numeric vector c(xmin, ymin, xmax, ymax) in degrees.
#' @param tile_dx,tile_dy Numeric; tile width/height in degrees.
#'
#' @return Data frame with columns x0, y0, x1, y1.
#'
#' @keywords internal
.make_tiles <- function(bbox_vec, tile_dx = 0.25, tile_dy = 0.25) {
  xmin <- bbox_vec[1]; ymin <- bbox_vec[2]
  xmax <- bbox_vec[3]; ymax <- bbox_vec[4]
  
  xs <- seq(xmin, xmax, by = tile_dx)
  ys <- seq(ymin, ymax, by = tile_dy)
  
  if (tail(xs, 1) < xmax) xs <- c(xs, xmax)
  if (tail(ys, 1) < ymax) ys <- c(ys, ymax)
  
  tiles <- expand.grid(
    x0 = head(xs, -1),
    y0 = head(ys, -1)
  )
  tiles$x1 <- pmin(tiles$x0 + tile_dx, xmax)
  tiles$y1 <- pmin(tiles$y0 + tile_dy, ymax)
  
  tiles
}


# =========================================================================
# Package setup and constants
# =========================================================================

required_pkgs <- c("keyring", "CDSE", "httr2", "terra", "sf", "move2", "dplyr")

missing_pkgs <- required_pkgs[!required_pkgs %in% installed.packages()[, "Package"]]
if (length(missing_pkgs) > 0) {
  message("Installing missing packages: ", paste(missing_pkgs, collapse = ", "))
  install.packages(missing_pkgs)
}

invisible(lapply(required_pkgs, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required but could not be loaded.")
  }
  library(pkg, character.only = TRUE)
}))

#' Sentinel Hub Process API evalscript for COPERNICUS_30 DEM
#'
#' @keywords internal
.dem_evalscript <- "
//VERSION=3
function setup() {
  return {
    input: [\"DEM\"],
    output: {
      id: \"default\",
      bands: 1,
      sampleType: SampleType.FLOAT32
    }
  };
}
function evaluatePixel(sample) {
  return [sample.DEM];
}
"


# =========================================================================
# Example usage
# =========================================================================

library(move2)
fishers <- mt_read(mt_example())

# Default 1 km buffer (safe for most cases)
res <- get_COP30_move2(fishers, 
                       tiling_strategy = "data-aware", 
                       key_service = "SentinelHUB")

# Increase buffer if you still get NAs near edges
res <- get_COP30_move2(fishers, 
                       tiling_strategy = "data-aware",
                       buffer_deg = 0.03,  # ~3 km buffer
                       key_service = "SentinelHUB")

fishers_annotated <- res$move2
head(fishers_annotated$dem_copernicus30)

# Inspect result:
fishers_with_dem <- res$move2
if (res$tiled) {
  cat("Used", length(res$dem_list), "tiles\n")
} else {
  cat("Single request\n")
}

leo <- readRDS("/home/kami/Documents/Teaching/Animove/2025CR/data/Leo-65545_seasons.rds")
leo
# Default 1 km buffer (safe for most cases)
res_leo <- get_COP30_move2(leo, 
                       tiling_strategy = "data-aware", 
                       key_service = "SentinelHUB")
leo_annotated <- res_leo$move2
head(leo_annotated$dem_copernicus30)

# Inspect result:
leo_with_dem <- res_leo$move2
if (res_leo$tiled) {
  cat("Used", length(res_leo$dem_list), "tiles\n")
} else {
  cat("Single request\n")
}
