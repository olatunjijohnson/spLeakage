# deleak_estimate(): turn the optimism gap into a bias-corrected ("de-leaked")
# accuracy estimate with a confidence interval -- the quantity practitioners report.
# Builds on the design-matched controlled scheme (the empirical optimism route) and
# adds a cluster (fold) bootstrap for uncertainty, plus the ability to correct a
# *reported* metric value without refitting (the meta-audit enabler, task #13).
# Grounded in docs/THEORY.md: the controlled error targets E_true (deployment error).

# Per-point loss contributions per fold (so the cluster bootstrap can resample them).
.fold_errors <- function(df, folds, predict_fun, response, metric) {
  y <- df[[response]]; loss <- .metric_spec(metric)$loss
  out <- lapply(folds, function(f) {
    if (!length(f$test) || !length(f$train)) return(numeric(0))
    loss(y[f$test], predict_fun(df[f$train, , drop = FALSE], df[f$test, , drop = FALSE]))
  })
  out[lengths(out) > 0L]
}

.agg_err <- function(err_list, metric) {
  .metric_spec(metric)$agg(unlist(err_list, use.names = FALSE))
}

# Cluster bootstrap over folds: resample whole folds with replacement.
.boot_metric <- function(err_list, metric, B) {
  K <- length(err_list)
  vapply(seq_len(B), function(b) .agg_err(err_list[sample.int(K, K, replace = TRUE)], metric),
         numeric(1))
}

# Target-matched control: among candidate schemes (random, block, and block buffered
# by increasing distances), pick the one whose SLI_rho against the declared target is
# closest to zero -- the scheme that best imitates deployment (NNDM matching as a
# selection criterion). Its error is the target-aware de-leaked estimate.
.matched_control <- function(data, target, response, coords, k, xy, info, dependence) {
  n <- nrow(xy$coords)
  if (is.null(dependence)) dependence <- estimate_dependence(data, response = response, coords = coords)
  tgt <- .as_prediction_target(target, data = data, coords = coords)
  diam <- sqrt(sum(apply(xy$coords, 2L, function(z) diff(range(z)))^2))
  set.seed(1)
  base_block <- .parse_split(spatial_block_cv(xy$coords, k), n)
  cand <- list(random = .parse_split(sample(rep_len(seq_len(k), n)), n), block = base_block)
  for (frac in c(0.05, 0.1, 0.15, 0.2, 0.25)) {
    cand[[sprintf("buffer%.2f", frac)]] <- .buffer_folds(xy$coords, base_block, frac * diam, info)
  }
  slis <- vapply(cand, function(f) tryCatch(
    detect_leakage(data, f, tgt, dependence = dependence, coords = coords)$SLI_rho,
    error = function(e) NA_real_), numeric(1))
  ok <- which(is.finite(slis))
  best <- ok[which.min(abs(slis[ok]))]
  list(folds = cand[[best]], sli = unname(slis[best]), label = names(cand)[best])
}

#' De-leaked (bias-corrected) accuracy estimate
#'
#' Reports the accuracy a model *would* achieve under a leakage-controlled,
#' design-matched validation -- the de-leaked estimate practitioners should report --
#' with a fold (cluster) bootstrap confidence interval. Optionally corrects a
#' user-supplied `reported` metric value (e.g. a published RMSE) using the estimated
#' optimism ratio, without refitting the original model.
#'
#' @inheritParams estimate_optimism
#' @param reported Optional reported metric value to correct (e.g. a published RMSE
#'   from a model not available for refitting). The de-leaked value is
#'   `reported * (controlled / user-CV)` error ratio.
#' @param target Optional [prediction_target()]. When supplied, the controlled
#'   scheme is **target-matched** (the candidate scheme whose `SLI_rho` against the
#'   target is closest to zero) instead of a fixed spatial block -- recommended when
#'   deployment is extrapolation, where a fixed block under-corrects.
#' @param n_boot Number of fold-bootstrap resamples for the interval.
#' @param level Confidence level for the interval.
#' @return An object of class `deleak_estimate`.
#' @examples
#' set.seed(1)
#' d <- data.frame(x = runif(150), y = runif(150), z = rnorm(150))
#' deleak_estimate(d, split = sample(rep_len(1:10, 150)), response = "z",
#'                 coords = c("x", "y"), n_boot = 100)
#' @export
deleak_estimate <- function(data, split, response, reported = NULL, target = NULL,
                            predict_fun = NULL, metric = c("rmse", "mae", "brier", "logloss", "poisson"),
                            control = c("block", "buffer"), coords = NULL,
                            dependence = NULL, k_control = NULL, buffer = NULL,
                            n_boot = 500L, level = 0.90) {
  metric <- match.arg(metric); control <- match.arg(control)
  xy <- .extract_coords(data, coords); n <- nrow(xy$coords)
  info <- list(crs = xy$crs, geographic = xy$geographic)
  base <- if (inherits(data, "sf")) sf::st_drop_geometry(data) else as.data.frame(data)
  df <- base; df$.x <- xy$coords[, 1]; df$.y <- xy$coords[, 2]
  if (!response %in% names(df)) stop(sprintf("Response column '%s' not found.", response))
  if (is.null(predict_fun)) predict_fun <- .idw_predictor(c(".x", ".y"), response)

  folds_user <- .parse_split(split, n)
  k <- if (!is.null(k_control)) k_control else max(length(folds_user), 5L)
  control_label <- control; matched_sli <- NA_real_
  if (!is.null(target)) {
    mc <- .matched_control(data, target, response, coords, k, xy, info, dependence)
    ctrl_folds <- mc$folds
    control_label <- paste0("target-matched (", mc$label, ")"); matched_sli <- mc$sli
  } else {
    ctrl_folds <- switch(
      control,
      block = .parse_split(spatial_block_cv(xy$coords, k), n),
      buffer = {
        buf <- buffer
        if (is.null(buf)) {
          if (is.null(dependence)) dependence <- estimate_dependence(data, response = response, coords = coords)
          buf <- dependence$practical_range
        }
        diameter <- sqrt(sum(apply(xy$coords, 2L, function(z) diff(range(z)))^2))
        if (buf > 0.25 * diameter) buf <- 0.25 * diameter
        .buffer_folds(xy$coords, .parse_split(spatial_block_cv(xy$coords, k), n), buf, info)
      })
  }

  eu <- .fold_errors(df, folds_user, predict_fun, response, metric)
  ec <- .fold_errors(df, ctrl_folds, predict_fun, response, metric)
  E_cv <- .agg_err(eu, metric); E_dl <- .agg_err(ec, metric)
  ratio <- if (E_cv > 0) E_dl / E_cv else NA_real_

  a <- (1 - level) / 2
  bc <- .boot_metric(ec, metric, n_boot)
  ratio_b <- bc / .boot_metric(eu, metric, n_boot)
  dl_ci <- unname(stats::quantile(bc, c(a, 1 - a)))
  ratio_ci <- unname(stats::quantile(ratio_b, c(a, 1 - a)))

  anchored <- anchored_ci <- NULL
  if (!is.null(reported)) { anchored <- reported * ratio; anchored_ci <- reported * ratio_ci }

  structure(
    list(metric = metric, control = control_label, matched_sli = matched_sli,
         reported_cv = E_cv, deleaked = E_dl, deleaked_ci = dl_ci,
         optimism = E_dl - E_cv, optimism_rel = if (E_dl > 0) (E_dl - E_cv) / E_dl else NA_real_,
         ratio = ratio, ratio_ci = ratio_ci,
         reported = reported, anchored = anchored, anchored_ci = anchored_ci,
         level = level, n = n),
    class = "deleak_estimate")
}

#' @export
print.deleak_estimate <- function(x, ...) {
  pct <- round(100 * (x$level), 0)
  cat("<deleak_estimate>\n")
  cat(sprintf("  metric / control     : %s / %s\n", toupper(x$metric), x$control))
  if (is.finite(x$matched_sli)) {
    cat(sprintf("  control match (SLI)  : %+.3f  (-> 0 is well matched to the target)\n", x$matched_sli))
  }
  cat(sprintf("  reported (your CV)   : %.4g\n", x$reported_cv))
  cat(sprintf("  de-leaked estimate   : %.4g   (%d%% CI %.4g, %.4g)\n",
              x$deleaked, pct, x$deleaked_ci[1], x$deleaked_ci[2]))
  cat(sprintf("  => reported accuracy inflated by %+.0f%% (error x%.2f)\n",
              100 * x$optimism_rel, x$ratio))
  if (!is.null(x$reported)) {
    cat(sprintf("  your reported %.4g -> de-leaked %.4g (%d%% CI %.4g, %.4g)\n",
                x$reported, x$anchored, pct, x$anchored_ci[1], x$anchored_ci[2]))
  }
  invisible(x)
}
