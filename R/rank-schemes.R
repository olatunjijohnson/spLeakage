# Data-driven scheme selection (contribution C3, upgraded). Instead of a static
# expert table, evaluate candidate CV schemes against the declared deployment target
# and rank them by how well their test->train geometry matches the deployment
# prediction->sample geometry (smallest |SLI_rho| = best matched). This operationalises
# the NNDM matching idea as a *selection criterion across schemes*, and makes the
# recommendation a computation rather than a lookup.

#' Rank candidate cross-validation schemes by deployment match
#'
#' Builds each candidate CV scheme, measures its Spatial Leakage Index against the
#' declared prediction target, and ranks them: the scheme whose `SLI_rho` is closest
#' to zero best imitates deployment (neither optimistic nor pessimistic).
#'
#' @param data An `sf` object, numeric matrix, or `data.frame`.
#' @param target A [prediction_target()] (or coordinates coercible to one).
#' @param response Name of the numeric response column (for the dependence model).
#' @param coords For non-`sf` input, the coordinate column names/indices.
#' @param k Number of folds for each candidate scheme.
#' @param schemes Candidate schemes: any of `"random"`, `"block"`, `"buffered"`.
#' @param dependence Optional precomputed [estimate_dependence()] object.
#' @param seed RNG seed for the random scheme.
#' @return An object of class `scheme_ranking`.
#' @examples
#' set.seed(1)
#' d <- data.frame(x = runif(150), y = runif(150), z = rnorm(150))
#' grid <- as.matrix(expand.grid(x = seq(0, 1, .2), y = seq(0, 1, .2)))
#' rank_cv_schemes(d, prediction_target(grid = grid, type = "grid"),
#'                 response = "z", coords = c("x", "y"))
#' @export
rank_cv_schemes <- function(data, target, response, coords = NULL, k = 10L,
                            schemes = c("random", "block", "buffered"),
                            dependence = NULL, seed = 1L) {
  schemes <- match.arg(schemes, several.ok = TRUE)
  xy <- .extract_coords(data, coords); n <- nrow(xy$coords)
  info <- list(crs = xy$crs, geographic = xy$geographic)
  if (is.null(dependence)) dependence <- estimate_dependence(data, response = response, coords = coords)
  tgt <- .as_prediction_target(target, data = data, coords = coords)

  build <- function(s) {
    set.seed(seed)
    switch(
      s,
      random   = .parse_split(sample(rep_len(seq_len(k), n)), n),
      block    = .parse_split(spatial_block_cv(xy$coords, k), n),
      buffered = {
        diam <- sqrt(sum(apply(xy$coords, 2L, function(z) diff(range(z)))^2))
        buf <- min(dependence$practical_range, 0.25 * diam)
        .buffer_folds(xy$coords, .parse_split(spatial_block_cv(xy$coords, k), n), buf, info)
      })
  }
  rows <- lapply(schemes, function(s) {
    lk <- detect_leakage(data, build(s), tgt, dependence = dependence, coords = coords)
    data.frame(scheme = s, sli_rho = lk$SLI_rho, abs_sli = abs(lk$SLI_rho),
               W = lk$W, verdict = .leak_verdict(lk$SLI_rho, lk$delta),
               stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, rows)
  tab <- tab[order(tab$abs_sli), , drop = FALSE]
  rownames(tab) <- NULL
  structure(list(ranking = tab, best = tab$scheme[1], target = tgt$type),
            class = "scheme_ranking")
}

#' @export
print.scheme_ranking <- function(x, ...) {
  cat(sprintf("<scheme_ranking>  target = '%s'\n", x$target))
  cat("  ranked by deployment match (|SLI_rho| -> 0 is best):\n")
  r <- x$ranking
  for (i in seq_len(nrow(r))) {
    cat(sprintf("    %s%-9s SLI_rho %+.3f  W %.3g  [%s]\n",
                if (i == 1L) "* " else "  ", r$scheme[i], r$sli_rho[i], r$W[i], r$verdict[i]))
  }
  cat(sprintf("  recommended: %s\n", x$best))
  invisible(x)
}
