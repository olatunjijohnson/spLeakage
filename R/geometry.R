# Geometry engine: coordinate extraction, CRS handling, nearest-neighbour distances.
# Internal helpers (not exported). See docs/METHOD-SLI.md notation.

# Geometry/CRS info from an sf crs (or NULL).
.geom_info <- function(crs = NA) {
  geographic <- FALSE
  if (inherits(crs, "crs")) {
    lonlat <- sf::st_is_longlat(crs)
    geographic <- isTRUE(lonlat)
  }
  list(crs = crs, geographic = geographic)
}

# Extract an n x 2 coordinate matrix + geometry info from sf, matrix, or data.frame.
# `coords` names/indices are required for non-sf input.
.extract_coords <- function(x, coords = NULL) {
  if (inherits(x, "sf") || inherits(x, "sfc")) {
    xy <- sf::st_coordinates(x)
    if (ncol(xy) < 2) stop("Could not extract 2D coordinates from sf object.")
    info <- .geom_info(sf::st_crs(x))
    return(list(coords = unname(xy[, 1:2, drop = FALSE]),
                crs = info$crs, geographic = info$geographic))
  }
  if (is.matrix(x) && is.numeric(x) && ncol(x) >= 2) {
    return(list(coords = unname(x[, 1:2, drop = FALSE]),
                crs = NA, geographic = FALSE))
  }
  if (is.data.frame(x)) {
    if (is.null(coords)) {
      stop("`coords` (column names or indices) is required for data.frame input.")
    }
    m <- as.matrix(x[, coords, drop = FALSE])
    storage.mode(m) <- "double"
    return(list(coords = unname(m), crs = NA, geographic = FALSE))
  }
  stop("Unsupported spatial input: provide an sf object, numeric matrix, or data.frame + `coords`.")
}

# Pull a response vector from sf/data.frame by column name.
.get_response <- function(data, response) {
  if (is.null(response)) return(NULL)
  if (inherits(data, "sf")) data <- sf::st_drop_geometry(data)
  if (!response %in% names(data)) {
    stop(sprintf("Response column '%s' not found in data.", response))
  }
  out <- data[[response]]
  if (!is.numeric(out)) {
    stop(sprintf("Response column '%s' must be numeric for variogram estimation.", response))
  }
  out
}

# 1-NN Euclidean distance from each query row to the reference set.
.knn1_euclid <- function(query, ref) {
  if (nrow(ref) == 0L) return(rep(NA_real_, nrow(query)))
  if (requireNamespace("FNN", quietly = TRUE)) {
    return(FNN::get.knnx(ref, query, k = 1)$nn.dist[, 1])
  }
  # Base fallback (O(n*m)); fine for MVP-scale data.
  apply(query, 1L, function(p) sqrt(min(colSums((t(ref) - p)^2))))
}

# Nearest-neighbour distance from each query point to a reference set, honouring
# the geometry: geodesic for geographic CRS, Euclidean otherwise.
.nn_dist <- function(query, ref, info) {
  if (nrow(ref) == 0L) return(rep(NA_real_, nrow(query)))
  if (!isTRUE(info$geographic)) {
    return(.knn1_euclid(query, ref))
  }
  q <- sf::st_as_sf(as.data.frame(query), coords = c(1L, 2L), crs = info$crs)
  r <- sf::st_as_sf(as.data.frame(ref), coords = c(1L, 2L), crs = info$crs)
  # st_nearest_feature emits an informational message for lon/lat (planar nearest);
  # the returned distance below is still geodesic. Silence the chatter.
  idx <- suppressMessages(sf::st_nearest_feature(q, r))
  as.numeric(suppressMessages(sf::st_distance(q, r[idx, ], by_element = TRUE)))
}

# Retained spatial correlation at each query point from its k nearest training
# points (noisy-OR combination 1 - prod(1 - rho(d_j))). Reduces exactly to
# rho(nearest) when k = 1, so the single-neighbour SLI is the k = 1 special case.
# Honest finding (docs/THEORY-RESULTS.md): k > 1 improves the *magnitude* match to
# exact GP optimism but *reduces* the correlation (noisy-OR over-counts mutually
# correlated neighbours and saturates), so k = 1 is the recommended default -- it is
# already near-sufficient for predicting optimism (cor 0.98). For k > 1 with
# geographic CRS, distances use planar coordinates.
.retained_corr <- function(query, ref, info, rho, k = 1L) {
  if (nrow(ref) == 0L) return(rep(NA_real_, nrow(query)))
  if (k <= 1L || nrow(ref) == 1L) return(rho(.nn_dist(query, ref, info)))
  k <- min(k, nrow(ref))
  if (requireNamespace("FNN", quietly = TRUE)) {
    D <- FNN::get.knnx(ref, query, k = k)$nn.dist
  } else {
    D <- t(apply(query, 1L, function(p) sort(sqrt(colSums((t(ref) - p)^2)))[seq_len(k)]))
  }
  1 - apply(1 - rho(D), 1L, prod)
}

# Trapezoidal integration.
.trapz <- function(x, y) {
  n <- length(x)
  if (n < 2L) return(0)
  sum(diff(x) * (y[-n] + y[-1L]) / 2)
}
