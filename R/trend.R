# Trend strength (contribution from the meta-audit, task #24). The meta-audit showed
# optimism has two sources: a short-range AUTOCORRELATION channel (the SLI) and a
# large-scale TREND channel that the SLI misses. Trend-dominated fields (e.g.
# rainfall, elevation) are highly optimistic under extrapolation deployment because a
# trend-blind model (IDW, RF, GP-without-drift) cannot extrapolate the trend, yet the
# SLI is silent. Trend strength predicts that optimism (cor = 0.75 across the audit
# corpus; |SLI| + trend gives R^2 = 0.62 vs 0.15 for SLI alone).

#' Trend strength: large-scale structure that resists extrapolation
#'
#' The proportion of response variance explained by a low-order polynomial of the
#' coordinates -- a measure of large-scale spatial trend. High trend strength means a
#' trend-blind model will extrapolate poorly, a leakage channel the autocorrelation
#' [detect_leakage()] index does not capture.
#'
#' @param data An `sf` object, numeric matrix, or `data.frame`.
#' @param response Name of the numeric response column.
#' @param coords For non-`sf` input, the coordinate column names/indices.
#' @param degree Polynomial degree of the coordinate trend surface (default 2).
#' @return A single number in `[0, 1]`: the trend `R^2`.
#' @examples
#' set.seed(1)
#' d <- data.frame(x = runif(100), y = runif(100))
#' d$z <- 3 * d$x + rnorm(100, 0, 0.2)        # strong x-trend
#' trend_strength(d, "z", coords = c("x", "y"))
#' @export
trend_strength <- function(data, response, coords = NULL, degree = 2L) {
  xy <- .extract_coords(data, coords)
  z <- .get_response(data, response)
  if (is.null(z)) stop("`response` is required for trend_strength().")
  s <- as.data.frame(scale(xy$coords)); names(s) <- c("x", "y"); s$z <- z
  f <- if (degree <= 1L) z ~ x + y
       else stats::as.formula(sprintf("z ~ poly(x, %d) + poly(y, %d) + x:y", degree, degree))
  r2 <- suppressWarnings(summary(stats::lm(f, data = s))$r.squared)
  if (!is.finite(r2)) 0 else r2
}
