# Joint spatiotemporal leakage (contribution C4). In space-time data, a test
# observation leaks if a training observation is close in BOTH space and time. The
# separate spatial (SLI) and temporal (lookahead) channels can each look safe while
# the joint space-time proximity still leaks. We measure leakage in a scaled
# space-time metric -- each axis divided by its dependence range, so space and time
# are comparable -- and report the joint index alongside the space-only and
# time-only indices to show they differ. Default deployment is a one-step-ahead
# forecast (predict the sampled locations at a future time).

# Per-test-point nearest-training distance (Euclidean, any dimension) under a split.
.fold_nn <- function(coords_mat, folds, n) {
  g <- rep(NA_real_, n)
  for (f in folds) {
    if (!length(f$test) || !length(f$train)) next
    g[f$test] <- .knn1_euclid(coords_mat[f$test, , drop = FALSE], coords_mat[f$train, , drop = FALSE])
  }
  g
}

#' Detect joint spatiotemporal leakage
#'
#' Measures leakage in a scaled space-time metric (each axis divided by its
#' dependence range) for a split, relative to a deployment target (by default a
#' one-step-ahead forecast). Reports the joint space-time index together with the
#' space-only and time-only indices, since the joint measure can differ from either.
#'
#' @param data An `sf` object, numeric matrix, or `data.frame`.
#' @param split A split specification (see [detect_leakage()]).
#' @param time Name of the time column (numeric or `Date`).
#' @param coords For non-`sf` input, the coordinate column names/indices.
#' @param sp_range,t_range Spatial / temporal dependence ranges used to scale the two
#'   axes (defaults: 30% of the spatial extent / temporal span). Supply estimated
#'   ranges for a calibrated metric.
#' @param horizon Forecast horizon for the default deployment target (defaults to the
#'   median temporal spacing).
#' @return An object of class `st_leakage`.
#' @export
detect_st_leakage <- function(data, split, time, coords = NULL,
                              sp_range = NULL, t_range = NULL, horizon = NULL) {
  xy <- .extract_coords(data, coords); n <- nrow(xy$coords)
  if (isTRUE(xy$geographic)) warning("Space-time distances use planar coordinates; ",
                                     "project geographic data for accuracy.", call. = FALSE)
  tt <- .time_vector(data, time)
  diam <- sqrt(sum(apply(xy$coords, 2L, function(z) diff(range(z)))^2))
  tspan <- diff(range(tt))
  if (is.null(sp_range)) sp_range <- 0.3 * diam
  if (is.null(t_range)) t_range <- 0.3 * tspan
  if (is.null(horizon)) horizon <- stats::median(diff(sort(unique(tt))))

  sp <- xy$coords / sp_range                       # scaled space (2 cols)
  tm <- matrix(tt / t_range, ncol = 1L)            # scaled time
  st <- cbind(sp, tm)                              # scaled space-time
  folds <- .parse_split(split, n)

  g_s <- .fold_nn(sp, folds, n); g_t <- .fold_nn(tm, folds, n); g_st <- .fold_nn(st, folds, n)
  tested <- which(!is.na(g_st))

  # Default deployment: forecast the sampled locations at t_max + horizon.
  p_sp <- xy$coords / sp_range
  p_tm <- matrix(rep((max(tt) + horizon) / t_range, n), ncol = 1L)
  p_st <- cbind(p_sp, p_tm)
  f_s  <- .knn1_euclid(p_sp, sp); f_t <- .knn1_euclid(p_tm, tm); f_st <- .knn1_euclid(p_st, st)

  sli <- function(g, f) .sli_areas(g, f)$A         # signed: + => test closer than deployment
  structure(
    list(sli_space = sli(g_s[tested], f_s), sli_time = sli(g_t[tested], f_t),
         sli_st = sli(g_st[tested], f_st),
         med_st_test = stats::median(g_st[tested]), med_st_deploy = stats::median(f_st),
         sp_range = sp_range, t_range = t_range, horizon = horizon, n_test = length(tested)),
    class = "st_leakage")
}

#' @export
print.st_leakage <- function(x, ...) {
  vd <- function(v) if (v > 0.02) "optimistic leakage" else if (v < -0.02) "pessimistic" else "matched"
  cat("<st_leakage>  (one-step-ahead forecast deployment)\n")
  cat(sprintf("  space-only SLI : %+.3f  [%s]\n", x$sli_space, vd(x$sli_space)))
  cat(sprintf("  time-only  SLI : %+.3f  [%s]\n", x$sli_time, vd(x$sli_time)))
  cat(sprintf("  JOINT space-time SLI : %+.3f  [%s]\n", x$sli_st, vd(x$sli_st)))
  cat(sprintf("  median space-time NN: test %.3f vs deployment %.3f (in dependence units)\n",
              x$med_st_test, x$med_st_deploy))
  invisible(x)
}
