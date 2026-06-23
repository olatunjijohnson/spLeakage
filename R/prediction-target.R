# The prediction-target API: the declared deployment geometry that SLI is read
# against (docs/METHOD-SLI.md section 5). The target is an explicit input, so
# leakage is a property of (split, data, declaration) -- never the split alone.

#' Declare the prediction target (deployment geometry)
#'
#' The Spatial Leakage Index is defined relative to where the model will actually
#' be used. This constructor declares that deployment geometry.
#'
#' @param data Sample data (used to derive an interpolation grid, and for CRS).
#' @param grid,newdata Explicit prediction locations (`sf`, matrix, or `data.frame`)
#'   for `type = "grid"` / `"newdata"`.
#' @param type One of `"grid"` (wall-to-wall mapping), `"newdata"` (supplied
#'   locations), or `"interpolation"` (unsampled locations within the sampled
#'   domain, generated from `data`).
#' @param coords For non-`sf` input, the coordinate column names/indices.
#' @param n Number of points to generate for `type = "interpolation"`.
#' @return An object of class `prediction_target` holding the prediction coordinates
#'   and geometry info.
#' @export
prediction_target <- function(data = NULL, grid = NULL, newdata = NULL,
                              type = c("grid", "newdata", "interpolation"),
                              coords = NULL, n = 5000L) {
  type <- match.arg(type)
  src <- switch(type, grid = grid, newdata = newdata, interpolation = NULL)

  if (type %in% c("grid", "newdata")) {
    if (is.null(src)) stop(sprintf("`%s` must be supplied for type = '%s'.", type, type))
    xy <- .extract_coords(src, coords)
    pc <- xy$coords; crs <- xy$crs; geo <- xy$geographic
  } else {
    if (is.null(data)) stop("`data` must be supplied for type = 'interpolation'.")
    xy <- .extract_coords(data, coords)
    bb <- apply(xy$coords, 2L, range)
    pc <- cbind(stats::runif(n, bb[1, 1], bb[2, 1]),
                stats::runif(n, bb[1, 2], bb[2, 2]))
    crs <- xy$crs; geo <- xy$geographic
  }
  structure(
    list(coords = unname(pc), crs = crs, geographic = geo, type = type),
    class = "prediction_target")
}

# Coerce convenience inputs (sf / matrix) into a prediction_target.
.as_prediction_target <- function(x, data = NULL, coords = NULL) {
  if (inherits(x, "prediction_target")) return(x)
  prediction_target(grid = x, type = "grid", coords = coords)
}

#' @export
print.prediction_target <- function(x, ...) {
  cat(sprintf("<prediction_target> type = '%s'\n", x$type))
  cat(sprintf("  prediction points : %d\n", nrow(x$coords)))
  cat(sprintf("  geographic CRS    : %s\n", isTRUE(x$geographic)))
  invisible(x)
}
