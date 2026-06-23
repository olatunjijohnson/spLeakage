# Optimal validation design (contribution ⑥). Inverts the SLI: given a trained model
# and a deployment target, place a budget of independent validation points so the
# estimate of true map accuracy is (i) unbiased -- the validation->training distance
# distribution matches the deployment's (NNDM matching as a *design* criterion), and
# (ii) minimum-variance -- distance-stratified Neyman allocation, with per-stratum
# error variance from the GP theory L(d) = V[1 - w^2 rho(d)^2] (docs/THEORY.md).
# Unifies leakage diagnosis with spatial survey design.

# Pick `n` spatially spread points from candidate coordinates (k-means centroids ->
# nearest candidate).
.select_spread <- function(coords, n) {
  if (n >= nrow(coords)) return(seq_len(nrow(coords)))
  km <- suppressWarnings(stats::kmeans(coords, centers = n, nstart = 3L, iter.max = 30L))
  vapply(seq_len(n), function(k) {
    idx <- which(km$cluster == k)
    idx[which.min(colSums((t(coords[idx, , drop = FALSE]) - km$centers[k, ])^2))]
  }, integer(1))
}

#' Design an optimal independent validation sample
#'
#' Given training locations and a deployment target, choose where to place a budget
#' of independent validation points so the resulting estimate of true map accuracy is
#' unbiased (its distance-to-training distribution matches deployment) and
#' minimum-variance (distance-stratified Neyman allocation using the GP error model).
#'
#' @param data Training data (`sf`, matrix, or `data.frame`).
#' @param target A [prediction_target()] (the deployment region).
#' @param budget Number of validation points to place.
#' @param response Name of the response column (for the dependence/error model).
#' @param coords For non-`sf` input, the coordinate column names/indices.
#' @param candidates Optional pool of candidate validation locations (default: the
#'   target locations).
#' @param allocation `"optimal"` (Neyman, default), `"proportional"`, or `"random"`.
#' @param strata Number of distance-to-training strata.
#' @param dependence Optional [estimate_dependence()] object.
#' @param seed RNG seed.
#' @return An object of class `validation_design`.
#' @examples
#' set.seed(1)
#' d <- data.frame(x = runif(80), y = runif(80), z = rnorm(80))
#' grid <- as.matrix(expand.grid(x = seq(0, 1, .1), y = seq(0, 1, .1)))
#' design_validation(d, prediction_target(grid = grid, type = "grid"),
#'                   budget = 20, response = "z", coords = c("x", "y"))
#' @export
design_validation <- function(data, target, budget, response = NULL, coords = NULL,
                              candidates = NULL, allocation = c("optimal", "proportional", "random"),
                              strata = 6L, dependence = NULL, seed = 1L) {
  allocation <- match.arg(allocation)
  set.seed(seed)
  xy <- .extract_coords(data, coords); info <- list(crs = xy$crs, geographic = xy$geographic)
  tgt <- .as_prediction_target(target, data = data, coords = coords)
  cand <- if (is.null(candidates)) tgt$coords else .extract_coords(candidates, coords)$coords

  d_cand <- .nn_dist(cand, xy$coords, info)
  d_tgt  <- .nn_dist(tgt$coords, xy$coords, info)

  if (is.null(dependence) && !is.null(response)) dependence <- estimate_dependence(data, response = response, coords = coords)
  Lfun <- if (!is.null(dependence)) {
    w <- dependence$signal_prop; rng <- dependence$range
    function(dd) 1 - w^2 * exp(-2 * dd / rng)               # GP pointwise MSE / V
  } else function(dd) dd / max(d_tgt + 1e-9)                # fallback: error ~ distance

  br <- unique(stats::quantile(d_tgt, seq(0, 1, length.out = strata + 1L)))
  sh_tgt <- cut(d_tgt, br, include.lowest = TRUE)
  sh_cand <- cut(d_cand, br, include.lowest = TRUE)
  W <- as.numeric(table(sh_tgt)) / length(d_tgt)            # deployment stratum weights
  L_h <- as.numeric(tapply(d_tgt, sh_tgt, function(dd) mean(Lfun(dd))))
  L_h[!is.finite(L_h)] <- mean(L_h, na.rm = TRUE)

  # Within-stratum variance of the per-point squared error (Gaussian: Var(e^2)=2 L^2).
  Sh2 <- 2 * L_h^2; Lbar <- sum(W * L_h)
  V_srs  <- sum(W * Sh2) + sum(W * (L_h - Lbar)^2)   # simple random validation
  V_prop <- sum(W * Sh2)                              # stratified, proportional
  V_opt  <- sum(W * sqrt(Sh2))^2                      # stratified, Neyman-optimal

  nlev <- length(W)
  if (allocation == "random") {
    idx <- sample(nrow(cand), min(budget, nrow(cand)))
    weights <- rep(1 / length(idx), length(idx))
    V_chosen <- V_srs
  } else {
    a <- if (allocation == "optimal") W * sqrt(Sh2) else W   # Neyman vs proportional
    n_h <- pmax(0L, as.integer(round(budget * a / sum(a))))
    lev <- levels(sh_cand); idx <- integer(0); weights <- numeric(0)
    for (h in seq_len(nlev)) {
      if (n_h[h] == 0L) next
      pool <- which(sh_cand == lev[h]); if (!length(pool)) next
      sel <- .select_spread(cand[pool, , drop = FALSE], min(n_h[h], length(pool)))
      sidx <- pool[sel]; idx <- c(idx, sidx)
      weights <- c(weights, rep(W[h] / length(sidx), length(sidx)))   # Horvitz-Thompson
    }
    weights <- weights / sum(weights)
    V_chosen <- if (allocation == "optimal") V_opt else V_prop
  }
  se_reduction <- 1 - sqrt(V_chosen / V_srs)          # vs simple random validation

  structure(
    list(locations = cand[idx, , drop = FALSE], index = idx, weights = weights,
         budget = budget, allocation = allocation, n_selected = length(idx),
         strata = data.frame(weight = round(W, 3), mean_error = round(L_h, 3),
                             n = as.integer(table(factor(sh_cand[idx], levels = levels(sh_cand))))),
         se_reduction_vs_random = se_reduction),
    class = "validation_design")
}

#' @export
print.validation_design <- function(x, ...) {
  cat(sprintf("<validation_design>  budget = %d, allocation = %s, placed = %d\n",
              x$budget, x$allocation, x$n_selected))
  if (x$allocation != "random") {
    cat(sprintf("  estimate true accuracy as a weighted mean of validation errors (weights returned)\n"))
    cat(sprintf("  SE of the accuracy estimate: %.0f%% lower than simple random validation (same budget)\n",
                100 * x$se_reduction_vs_random))
  }
  cat("  distance strata (deployment weight x mean error x points placed):\n")
  print(x$strata, row.names = FALSE)
  invisible(x)
}
