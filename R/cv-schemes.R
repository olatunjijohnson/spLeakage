# Cross-validation scheme construction and a model-agnostic CV runner, used by
# estimate_optimism(). See docs/METHOD-EMULATOR.md section 1.

#' Spatial block cross-validation folds
#'
#' Assigns observations to `k` spatially contiguous folds by k-means clustering of
#' the coordinates. A leakage-controlled scheme for comparison against a user's
#' (often random) split.
#'
#' @param data An `sf` object, numeric coordinate matrix, or `data.frame`.
#' @param k Number of spatial folds.
#' @param coords For non-`sf` input, the coordinate column names/indices.
#' @return An integer fold-id vector of length `nrow(data)`.
#' @export
spatial_block_cv <- function(data, k = 10L, coords = NULL) {
  xy <- .extract_coords(data, coords)
  k <- min(k, nrow(unique(xy$coords)))
  km <- suppressWarnings(stats::kmeans(xy$coords, centers = k, nstart = 5L, iter.max = 50L))
  as.integer(km$cluster)
}

# Buffered folds: from a base fold structure, drop training points within `buffer`
# of any test point (buffered/dead-zone CV; Ploton et al. 2020).
.buffer_folds <- function(coords, folds, buffer, info) {
  lapply(folds, function(f) {
    if (!length(f$test) || !length(f$train)) return(f)
    d <- .nn_dist(coords[f$train, , drop = FALSE], coords[f$test, , drop = FALSE], info)
    list(test = f$test, train = f$train[d > buffer])
  })
}

# Model-agnostic out-of-fold prediction. `predict_fun(train_df, test_df)` returns
# numeric predictions for the test rows.
.cv_predict <- function(df, folds, predict_fun) {
  pred <- rep(NA_real_, nrow(df))
  for (f in folds) {
    if (!length(f$test) || !length(f$train)) next
    pred[f$test] <- predict_fun(df[f$train, , drop = FALSE], df[f$test, , drop = FALSE])
  }
  pred
}

# Inverse-distance-weighting predictor factory (purely spatial; the default learner
# so estimate_optimism works with no covariates and clearly exhibits leakage).
.idw_predictor <- function(coord_cols, response, power = 2) {
  function(train, test) {
    tr <- as.matrix(train[, coord_cols, drop = FALSE])
    te <- as.matrix(test[, coord_cols, drop = FALSE])
    y <- train[[response]]
    apply(te, 1L, function(p) {
      d <- sqrt(colSums((t(tr) - p)^2))
      hit <- which(d == 0)
      if (length(hit)) return(y[hit[1L]])
      w <- 1 / d^power
      sum(w * y) / sum(w)
    })
  }
}

# Error metrics, including proper scoring rules for non-Gaussian responses. Each
# spec gives a per-point `loss` (so the fold/cluster bootstrap can resample it) and
# an `agg` aggregator. All are oriented so lower = better.
#   rmse/mae  : continuous;  brier/logloss : binary (p in [0,1]);  poisson : counts.
.metric_spec <- function(metric) {
  switch(
    metric,
    rmse    = list(loss = function(o, p) (o - p)^2, agg = function(x) sqrt(mean(x))),
    mae     = list(loss = function(o, p) abs(o - p), agg = mean),
    brier   = list(loss = function(o, p) (o - p)^2, agg = mean),
    logloss = list(loss = function(o, p) {
      p <- pmin(pmax(p, 1e-6), 1 - 1e-6); -(o * log(p) + (1 - o) * log(1 - p))
    }, agg = mean),
    poisson = list(loss = function(o, p) {
      p <- pmax(p, 1e-9); 2 * (p - o + ifelse(o > 0, o * log(o / p), 0))
    }, agg = mean),
    stop(sprintf("Unknown metric '%s'.", metric)))
}

# Full metric as a single function(observed, predicted).
.metric_fun <- function(metric) {
  s <- .metric_spec(metric); function(o, p) s$agg(s$loss(o, p))
}
