# Feature-space leakage (contribution C4). Leakage is not only geographic: a test
# point that is close to training points in COVARIATE space is also evaluated more
# easily than a deployment point that lies in a sparsely-sampled region of feature
# space (cf. the area of applicability, Meyer & Pebesma 2021). This mirrors the
# geographic SLI but measures nearest-neighbour reach in standardised covariate space.

#' Detect feature-space (covariate) leakage
#'
#' Compares how close test points are to training points in standardised covariate
#' space against how close the deployment locations are to the sample. A positive
#' index means the test set is easier (closer in feature space) than deployment --
#' optimistic feature-space leakage / hidden extrapolation at deployment.
#'
#' @param data An `sf` object or `data.frame` containing the covariate columns.
#' @param split A split specification (see [detect_leakage()]).
#' @param covariates Character vector of covariate column names.
#' @param newdata Optional deployment data with the same covariate columns. Without
#'   it, only the test-set area-of-applicability is reported (no signed index).
#' @param coords Unused for the index (kept for API symmetry).
#' @return An object of class `feature_leakage`.
#' @export
detect_feature_leakage <- function(data, split, covariates, newdata = NULL,
                                   coords = NULL) {
  df <- if (inherits(data, "sf")) sf::st_drop_geometry(data) else as.data.frame(data)
  if (!all(covariates %in% names(df))) stop("Some `covariates` not found in data.")
  X <- as.matrix(df[, covariates, drop = FALSE]); storage.mode(X) <- "double"
  n <- nrow(X)
  center <- colMeans(X); scl <- apply(X, 2L, stats::sd); scl[scl == 0] <- 1
  Xs <- scale(X, center = center, scale = scl)

  folds <- .parse_split(split, n)
  gobs <- rep(NA_real_, n)
  for (k in seq_along(folds)) {
    f <- folds[[k]]
    if (!length(f$test) || !length(f$train)) next
    gobs[f$test] <- .knn1_euclid(Xs[f$test, , drop = FALSE], Xs[f$train, , drop = FALSE])
  }
  tested <- which(!is.na(gobs))

  # Feature-space scale: mean nearest-neighbour distance among the sample.
  if (requireNamespace("FNN", quietly = TRUE)) {
    selfnn <- FNN::get.knn(Xs, k = 1L)$nn.dist[, 1]
  } else {
    selfnn <- vapply(seq_len(n), function(i)
      min(sqrt(colSums((t(Xs[-i, , drop = FALSE]) - Xs[i, ])^2))), numeric(1))
  }
  avg <- mean(selfnn)
  thr <- as.numeric(stats::quantile(selfnn / avg, 0.95))
  di_test <- gobs[tested] / avg
  aoa_test <- mean(di_test <= thr)            # fraction of test points within applicability

  feature_sli <- NA_real_; mean_pred_di <- NA_real_
  if (!is.null(newdata)) {
    ndf <- if (inherits(newdata, "sf")) sf::st_drop_geometry(newdata) else as.data.frame(newdata)
    Xn <- scale(as.matrix(ndf[, covariates, drop = FALSE]), center = center, scale = scl)
    fpred <- .knn1_euclid(Xn, Xs)
    ar <- .sli_areas(gobs[tested], fpred)
    feature_sli <- ar$A / avg                 # signed; positive = optimistic feature leakage
    mean_pred_di <- mean(fpred / avg)
  }
  structure(
    list(feature_sli = feature_sli, n_test = length(tested),
         mean_test_di = mean(di_test), mean_pred_di = mean_pred_di,
         aoa_test = aoa_test, aoa_threshold = thr, avg_dist = avg,
         covariates = covariates, has_newdata = !is.null(newdata)),
    class = "feature_leakage")
}

#' @export
print.feature_leakage <- function(x, ...) {
  cat("<feature_leakage>\n")
  cat(sprintf("  covariates        : %s\n", paste(x$covariates, collapse = ", ")))
  cat(sprintf("  test in feature-AOA : %.0f%% (interpolation in covariate space)\n",
              100 * x$aoa_test))
  if (x$has_newdata) {
    dir <- if (x$feature_sli > 0.05) "OPTIMISTIC feature leakage"
           else if (x$feature_sli < -0.05) "deployment easier than test"
           else "matched"
    cat(sprintf("  feature SLI (signed): %+.3f  [%s]\n", x$feature_sli, dir))
    cat(sprintf("  mean NN reach (rel) : test %.2f vs deployment %.2f\n",
                x$mean_test_di, x$mean_pred_di))
  } else {
    cat("  (supply `newdata` for a signed feature-leakage index)\n")
  }
  invisible(x)
}
