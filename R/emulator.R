# The optimism emulator (contribution C2b): a fast, no-refit estimate of optimism
# from cheap data/split/target features, calibrated once on a simulation study with
# known ground-truth optimism. The emulator is conditioned on the model class and
# response type (categorical), with a continuous-feature k-NN + area-of-applicability
# guard so it refuses to extrapolate. See docs/METHOD-EMULATOR.md.

# Continuous (data-derived) features. Computed IDENTICALLY in the simulation study
# and at predict time -- this is the contract that keeps the emulator coherent.
# (sli_d is intentionally excluded: A/phi is numerically unstable when the fitted
# range is tiny, and it is uninformative relative to the dependence-form sli_rho.)
.EMU_FEATURES <- c("n", "range_rel", "signal", "sli_rho", "delta", "nn_index")

.emulator_features <- function(data, split, target, response, coords = NULL,
                               dependence = NULL) {
  if (is.null(dependence)) {
    dependence <- estimate_dependence(data, response = response, coords = coords)
  }
  lk <- detect_leakage(data, split, target, dependence = dependence, coords = coords)
  xy <- .extract_coords(data, coords)
  diam <- sqrt(sum(apply(xy$coords, 2L, function(z) diff(range(z)))^2))
  nnidx <- .nn_index(xy$coords, list(geographic = xy$geographic))
  c(n = lk$n,
    range_rel = dependence$practical_range / diam,
    signal = dependence$signal_prop,
    sli_rho = lk$SLI_rho,
    sli_d = lk$SLI_d,
    delta = lk$delta,
    nn_index = nnidx)
}

# Detect the response distribution from the values.
.detect_response_type <- function(y) {
  y <- y[is.finite(y)]
  if (all(y %in% c(0, 1))) "binomial"
  else if (all(y >= 0 & abs(y - round(y)) < 1e-8)) "poisson"
  else "gaussian"
}

# Categorical levels and the design-matrix encoding used by the gradient-boosted
# learner. Built identically in fit and predict (the contract).
.EMU_MODELS <- c("idw", "rf", "gam")
.EMU_RESP <- c("gaussian", "poisson", "binomial")
.emu_matrix <- function(df) {
  X <- as.matrix(df[, .EMU_FEATURES, drop = FALSE]); storage.mode(X) <- "double"
  for (m in .EMU_MODELS) X <- cbind(X, as.numeric(df$model == m))
  for (r in .EMU_RESP)   X <- cbind(X, as.numeric(df$response_type == r))
  colnames(X) <- c(.EMU_FEATURES, paste0("model_", .EMU_MODELS), paste0("resp_", .EMU_RESP))
  X
}

.emu_params <- list(objective = "reg:squarederror", eta = 0.05, max_depth = 4L,
                    subsample = 0.8, colsample_bytree = 0.8, min_child_weight = 5)
.EMU_NROUNDS <- 250L

# Fit the emulator from a feature/label table `tab` (.EMU_FEATURES + `model`,
# `response_type`, `optimism_rel`). Point estimate = gradient-boosted trees
# (xgboost); prediction intervals = split-conformal (5-fold out-of-fold residuals,
# guaranteeing ~90% marginal coverage); area-of-applicability = instance-based
# distance in standardised continuous feature space. Used by data-raw.
.fit_optimism_emulator <- function(tab, alpha = 0.10) {
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("Building the emulator requires the 'xgboost' package.", call. = FALSE)
  }
  ok <- stats::complete.cases(tab[, .EMU_FEATURES]) & is.finite(tab$optimism_rel)
  tab <- tab[ok, , drop = FALSE]
  if (nrow(tab) < 50L) stop("Too few simulation rows to fit the emulator.")
  X <- .emu_matrix(tab); y <- tab$optimism_rel

  # 5-fold out-of-fold residuals, grouped by config (a config's correlated rows held
  # out together), so coverage is nominal for unseen *datasets*, not unseen rows.
  set.seed(123)
  gid <- if (!is.null(tab$config)) tab$config else seq_len(nrow(X))
  ug <- unique(gid); fold <- sample(rep_len(1:5, length(ug)))[match(gid, ug)]
  oof <- numeric(nrow(X))
  for (k in 1:5) {
    tr <- fold != k
    mk <- xgboost::xgb.train(.emu_params, xgboost::xgb.DMatrix(X[tr, ], label = y[tr]),
                             nrounds = .EMU_NROUNDS, verbose = 0)
    oof[!tr] <- stats::predict(mk, X[!tr, ])
  }
  resid <- abs(y - oof)

  m <- xgboost::xgb.train(.emu_params, xgboost::xgb.DMatrix(X, label = y),
                          nrounds = .EMU_NROUNDS, verbose = 0)
  model_raw <- xgboost::xgb.save.raw(m)

  # Area of applicability + a LOCAL difficulty scale sigma(x) = mean residual of the
  # k_sigma nearest neighbours in standardised feature space (category-conditional).
  k_sigma <- 25L
  Xc <- as.matrix(tab[, .EMU_FEATURES])
  center <- colMeans(Xc); scl <- apply(Xc, 2L, stats::sd); scl[scl == 0] <- 1
  Xs <- scale(Xc, center = center, scale = scl)
  avg_dist <- mean(stats::dist(Xs))
  grp <- interaction(tab$model, tab$response_type, drop = TRUE)
  di <- rep(NA_real_, nrow(Xs)); sigma <- rep(NA_real_, nrow(Xs))
  for (g in levels(grp)) {
    idx <- which(grp == g); if (length(idx) < 2L) next
    D <- as.matrix(stats::dist(Xs[idx, , drop = FALSE])); diag(D) <- Inf
    di[idx] <- apply(D, 1L, min) / avg_dist
    ks <- min(k_sigma, length(idx) - 1L)
    sigma[idx] <- vapply(seq_along(idx), function(j) mean(resid[idx[order(D[j, ])[seq_len(ks)]]]), numeric(1))
  }
  sigma_floor <- 0.25 * stats::median(resid)
  sigma[!is.finite(sigma)] <- stats::median(resid)
  sigma <- pmax(sigma, sigma_floor)
  # Normalised-conformal multiplier: interval = pred +/- conf_q * sigma(x_query).
  conf_q <- as.numeric(stats::quantile(resid / sigma, 1 - alpha))
  qd <- stats::quantile(di[is.finite(di)], c(0.25, 0.75))
  thr <- as.numeric(qd[2] + 1.5 * (qd[2] - qd[1]))

  structure(
    list(model_raw = model_raw, conf_q = conf_q, alpha = alpha,
         resid = resid, k_sigma = k_sigma, sigma_floor = sigma_floor,
         Xs = unclass(Xs), model = as.character(tab$model),
         response_type = as.character(tab$response_type),
         center = center, scale = scl, features = .EMU_FEATURES,
         avg_dist = avg_dist, aoa_threshold = thr, n_train = nrow(Xs),
         learner = "xgboost + normalized conformal"),
    class = "optimism_emulator")
}

.get_emulator <- function() {
  e <- get0("optimism_emulator", envir = asNamespace("spLeakage"), inherits = FALSE)
  if (is.null(e)) {
    stop("Optimism emulator not built yet. Run data-raw/simulate_optimism.R to ",
         "generate it, or use estimate_optimism() (the empirical route).", call. = FALSE)
  }
  e
}

# Restore the xgboost booster from its raw serialisation.
.emu_booster <- function(emulator) xgboost::xgb.load.raw(emulator$model_raw)

# Core predictor. `booster` may be passed to avoid repeated deserialisation.
.emu_predict_one <- function(emulator, x_cont, model, response_type, booster = NULL) {
  x_std <- (x_cont[emulator$features] - emulator$center) / emulator$scale
  match_full <- emulator$model == model & emulator$response_type == response_type
  level <- "model+response"; idx <- which(match_full)
  if (length(idx) < 10L) { idx <- which(emulator$model == model); level <- "model" }
  if (length(idx) < 10L) { idx <- seq_along(emulator$model); level <- "none" }
  d <- sqrt(colSums((t(emulator$Xs[idx, , drop = FALSE]) - x_std)^2))
  di <- min(d) / emulator$avg_dist
  in_aoa <- (di <= emulator$aoa_threshold) && level != "none"

  if (is.null(booster)) booster <- .emu_booster(emulator)
  qdf <- data.frame(as.list(x_cont[emulator$features]),
                    model = model, response_type = response_type,
                    stringsAsFactors = FALSE)
  pred <- as.numeric(stats::predict(booster, .emu_matrix(qdf)))
  # Local difficulty scale: mean residual of the k_sigma nearest training points.
  ks <- min(emulator$k_sigma %||% 25L, length(d))
  sigma <- mean(emulator$resid[idx[order(d)[seq_len(ks)]]])
  sigma <- max(sigma, emulator$sigma_floor %||% 0)
  hw <- emulator$conf_q * sigma                     # normalized-conformal half-width
  ci <- c(pred - hw, pred + hw)
  if (!in_aoa) { pred <- NA_real_; ci <- c(NA_real_, NA_real_) }
  list(optimism_rel = pred, optimism_ci = ci, in_aoa = in_aoa, di = di, match = level)
}

#' Predict optimism from the calibrated emulator (fast, no model refitting)
#'
#' Estimates the optimism of a split using the pre-trained emulator and cheap
#' data/split/target features -- the fast counterpart to [estimate_optimism()]. The
#' estimate is conditioned on the model class and response type. If the query falls
#' outside the emulator's area of applicability (in feature or category space) the
#' function refuses a point estimate (`NA`) rather than extrapolating.
#'
#' @inheritParams detect_leakage
#' @param model The learner whose optimism to estimate: `"idw"` (default), `"rf"`,
#'   or `"gam"`.
#' @param response_type `"auto"` (detect from the response, default), `"gaussian"`,
#'   `"poisson"`, or `"binomial"`.
#' @param emulator An `optimism_emulator`; defaults to the one shipped with the package.
#' @return An object of class `optimism_prediction`.
#' @export
predict_optimism <- function(data, split, target, response = NULL,
                             model = c("idw", "rf", "gam"),
                             response_type = c("auto", "gaussian", "poisson", "binomial"),
                             dependence = NULL, coords = NULL, emulator = NULL) {
  model <- match.arg(model); response_type <- match.arg(response_type)
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("predict_optimism() needs the 'xgboost' package (the emulator's learner). ",
         "Install it, or use estimate_optimism() (the empirical route).", call. = FALSE)
  }
  if (is.null(emulator)) emulator <- .get_emulator()
  if (response_type == "auto") {
    y <- .get_response(data, response)
    response_type <- if (is.null(y)) "gaussian" else .detect_response_type(y)
  }
  feats <- .emulator_features(data, split, target, response = response,
                              coords = coords, dependence = dependence)
  res <- .emu_predict_one(emulator, feats, model, response_type)
  if (!res$in_aoa) {
    warning("Query is outside the emulator's area of applicability; refusing a ",
            "point estimate. Use estimate_optimism() instead.", call. = FALSE)
  }
  structure(
    c(res, list(aoa_threshold = emulator$aoa_threshold, model = model,
                response_type = response_type, features = feats)),
    class = "optimism_prediction")
}

#' @export
print.optimism_emulator <- function(x, ...) {
  cat("<optimism_emulator>\n")
  cat(sprintf("  trained on    : %d simulation rows\n", x$n_train))
  cat(sprintf("  models        : %s\n", paste(sort(unique(x$model)), collapse = ", ")))
  cat(sprintf("  responses     : %s\n", paste(sort(unique(x$response_type)), collapse = ", ")))
  cat(sprintf("  features      : %s\n", paste(x$features, collapse = ", ")))
  cat(sprintf("  learner       : %s\n", x$learner %||% "instance-based"))
  cat(sprintf("  AOA threshold : %.3f   conformal half-width: %.3f\n",
              x$aoa_threshold, x$conf_q %||% NA_real_))
  v <- attr(x, "validation")
  if (!is.null(v)) {
    cat(sprintf("  held-out val. : R2=%.2f, RMSE=%.3f, in-AOA=%.0f%%, 90%%-cover=%.0f%%\n",
                v$R2, v$rmse, 100 * v$aoa_coverage, 100 * v$interval_coverage))
  }
  invisible(x)
}

#' @export
print.optimism_prediction <- function(x, ...) {
  cat("<optimism_prediction> (emulator)\n")
  cat(sprintf("  model / response : %s / %s  [match: %s]\n", x$model, x$response_type, x$match))
  if (is.na(x$optimism_rel)) {
    cat(sprintf("  optimism         : NA  [outside AOA: DI %.2f > %.2f]\n",
                x$di, x$aoa_threshold))
  } else {
    cat(sprintf("  optimism         : %+.1f%%  (90%% CI %+.1f%%, %+.1f%%)  [in AOA]\n",
                100 * x$optimism_rel, 100 * x$optimism_ci[1], 100 * x$optimism_ci[2]))
  }
  invisible(x)
}
