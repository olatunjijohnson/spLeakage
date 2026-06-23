# Temporal / spatiotemporal leakage (contribution C4). A random split of time-stamped
# data trains on the future to predict the past ("lookahead" leakage), inflating
# performance relative to the real forecasting task. This audits a split for
# lookahead and time proximity; temporal_kfold() provides the forward-chaining fix.

.time_vector <- function(data, time) {
  df <- if (inherits(data, "sf")) sf::st_drop_geometry(data) else as.data.frame(data)
  if (!time %in% names(df)) stop(sprintf("Time column '%s' not found.", time))
  as.numeric(df[[time]])
}

#' Detect temporal (lookahead) leakage in a split
#'
#' Flags test observations whose fold's training set contains *later* time points
#' (the model is trained on the future to predict the past) and summarises how close
#' in time test points are to their training data.
#'
#' @param data An `sf` object or `data.frame`.
#' @param split A split specification (see [detect_leakage()]).
#' @param time Name of the time column (numeric or `Date`).
#' @param coords Unused (kept for API symmetry).
#' @return An object of class `temporal_leakage`.
#' @export
detect_temporal_leakage <- function(data, split, time, coords = NULL) {
  t <- .time_vector(data, time)
  n <- length(t)
  folds <- .parse_split(split, n)
  lookahead <- logical(n); tnn <- rep(NA_real_, n); tested <- logical(n)
  for (f in folds) {
    if (!length(f$test) || !length(f$train)) next
    tr <- t[f$train]
    for (i in f$test) {
      lookahead[i] <- any(tr > t[i])
      tnn[i] <- min(abs(tr - t[i]))
    }
    tested[f$test] <- TRUE
  }
  idx <- which(tested)
  structure(
    list(lookahead_frac = mean(lookahead[idx]), n_test = length(idx),
         n_lookahead = sum(lookahead[idx]), mean_time_nn = mean(tnn[idx]),
         time = time),
    class = "temporal_leakage")
}

#' Forward-chaining temporal cross-validation (the fix for lookahead leakage)
#'
#' Splits the data into `k` time-ordered blocks and builds folds where each block is
#' tested using only earlier blocks for training -- no future information leaks into
#' training. Returns a list of folds usable directly by [detect_leakage()] etc.
#'
#' @inheritParams detect_temporal_leakage
#' @param k Number of time blocks.
#' @return A list of folds, each `list(test=, train=)`.
#' @export
temporal_kfold <- function(data, k = 5L, time, coords = NULL) {
  t <- .time_vector(data, time)
  n <- length(t)
  ord <- order(t)
  block <- as.integer(cut(seq_len(n), breaks = k, labels = FALSE))
  folds <- list()
  for (i in 2:k) {
    test <- ord[block == i]; train <- ord[block < i]
    if (length(test) && length(train)) {
      folds[[length(folds) + 1L]] <- list(test = test, train = train)
    }
  }
  folds
}

#' @export
print.temporal_leakage <- function(x, ...) {
  cat("<temporal_leakage>\n")
  cat(sprintf("  time column       : %s   |  n test = %d\n", x$time, x$n_test))
  cat(sprintf("  lookahead leakage : %d / %d (%.1f%%) test points trained on the future\n",
              x$n_lookahead, x$n_test, 100 * x$lookahead_frac))
  cat(sprintf("  mean time-to-nearest-training : %.4g\n", x$mean_time_nn))
  if (x$lookahead_frac > 0) cat("  [!] fix: temporal_kfold() (forward-chaining CV)\n")
  invisible(x)
}
