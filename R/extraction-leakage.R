# Covariate-extraction-overlap leakage (contribution C4). Spatial covariates are
# often built from a neighbourhood window around each location -- focal/moving-window
# statistics, buffer summaries, kernel densities, distance-to-feature. If the
# extraction window of a TEST point overlaps that of a TRAINING point, the test
# covariate value encodes data from the training neighbourhood: leakage through the
# covariate *construction*, even when the points themselves are far enough apart for
# a spatial CV fold. This is distinct from autocorrelation leakage (the SLI) and from
# feature-space leakage (covariate *values*); it is about shared raw data in the
# extraction. Common and under-recognised in remote sensing / species distribution
# modelling. The fix is a train/test buffer of at least twice the extraction radius.

#' Detect covariate-extraction-overlap leakage
#'
#' Flags test points whose covariate-extraction window overlaps a training point's,
#' given the extraction radius (focal-window radius, buffer distance, or kernel
#' bandwidth) used to build the spatial covariates.
#'
#' @param data An `sf` object, numeric matrix, or `data.frame`.
#' @param split A split specification (see [detect_leakage()]).
#' @param radius The extraction-window radius / buffer / bandwidth, in coordinate
#'   units (geodesic metres for geographic CRS).
#' @param coords For non-`sf` input, the coordinate column names/indices.
#' @return An object of class `extraction_leakage`.
#' @examples
#' set.seed(1)
#' d <- data.frame(x = runif(100), y = runif(100), z = rnorm(100))
#' # covariates were extracted with a 0.1-unit focal window:
#' detect_extraction_leakage(d, split = sample(rep_len(1:5, 100)),
#'                           radius = 0.1, coords = c("x", "y"))
#' @export
detect_extraction_leakage <- function(data, split, radius, coords = NULL) {
  if (!is.numeric(radius) || radius <= 0) stop("`radius` must be a positive number.")
  xy <- .extract_coords(data, coords); n <- nrow(xy$coords)
  info <- list(crs = xy$crs, geographic = xy$geographic)
  folds <- .parse_split(split, n)
  gobs <- rep(NA_real_, n)
  for (f in folds) {
    if (!length(f$test) || !length(f$train)) next
    gobs[f$test] <- .nn_dist(xy$coords[f$test, , drop = FALSE],
                             xy$coords[f$train, , drop = FALSE], info)
  }
  tested <- which(!is.na(gobs)); g <- gobs[tested]
  structure(
    list(radius = radius, n_test = length(tested),
         n_overlap = sum(g < 2 * radius), frac_overlap = mean(g < 2 * radius),
         n_contains = sum(g < radius), frac_contains = mean(g < radius),
         recommended_buffer = 2 * radius),
    class = "extraction_leakage")
}

#' @export
print.extraction_leakage <- function(x, ...) {
  cat(sprintf("<extraction_leakage>  extraction radius = %.4g\n", x$radius))
  cat(sprintf("  test windows overlapping a training window (< 2r): %d / %d (%.1f%%)\n",
              x$n_overlap, x$n_test, 100 * x$frac_overlap))
  cat(sprintf("  test windows containing a training location (< r): %d / %d (%.1f%%)\n",
              x$n_contains, x$n_test, 100 * x$frac_contains))
  if (x$frac_overlap > 0) {
    cat(sprintf("  [!] fix: separate train/test by a buffer >= %.4g (= 2 x radius)\n",
                x$recommended_buffer))
  }
  invisible(x)
}
