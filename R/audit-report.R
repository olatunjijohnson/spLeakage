# audit_workflow() and report_leakage(): assemble the diagnostics into a
# submission-ready scorecard. See docs/VISION.md C5 and G.

# Letter grade from the (absolute) leakage index.
.grade <- function(sli) {
  a <- abs(sli)
  if (a < 0.05) "A" else if (a < 0.10) "B" else if (a < 0.20) "C" else
    if (a < 0.35) "D" else "F"
}

#' Audit a cross-validation workflow for spatial leakage
#'
#' Runs [detect_leakage()] and adds data-hygiene checks (duplicated coordinates,
#' documented CRS) to produce a scorecard.
#'
#' @inheritParams detect_leakage
#' @param group Optional grouping column for [detect_group_leakage()] (defaults to
#'   grouping by coordinates).
#' @return An object of class `workflow_audit`.
#' @export
audit_workflow <- function(data, split, target, response = NULL,
                           dependence = NULL, coords = NULL, group = NULL) {
  lk <- detect_leakage(data, split, target, response = response,
                       dependence = dependence, coords = coords)
  gl <- detect_group_leakage(data, split, group = group, coords = coords)
  xy <- .extract_coords(data, coords)
  trend <- if (!is.null(response)) tryCatch(trend_strength(data, response, coords),
                                            error = function(e) NA_real_) else NA_real_
  flags <- list(
    leakage_grade = .grade(lk$SLI_rho),
    sli_rho = lk$SLI_rho,
    trend_strength = trend,
    duplicated_coords = sum(duplicated(xy$coords)),
    group_leak_frac = gl$frac_leaked,
    crs_documented = inherits(xy$crs, "crs") && !is.na(xy$crs))
  structure(list(leakage = lk, group_leakage = gl, flags = flags),
            class = "workflow_audit")
}

#' @export
print.workflow_audit <- function(x, ...) {
  f <- x$flags
  cat("<workflow_audit>\n")
  cat(sprintf("  leakage grade     : %s  (SLI_rho = %+.3f)\n", f$leakage_grade, f$sli_rho))
  cat(sprintf("  %s autocorrelation leakage (random-split optimism)\n",
              if (f$sli_rho > 0.1) "[!]" else "[ok]"))
  if (is.finite(f$trend_strength)) {
    cat(sprintf("  %s trend strength: %.2f (extrapolation-optimism risk a trend-blind model misses)\n",
                if (f$trend_strength > 0.5) "[!]" else "[ok]", f$trend_strength))
  }
  cat(sprintf("  %s duplicated coordinates: %d\n",
              if (f$duplicated_coords > 0) "[!]" else "[ok]", f$duplicated_coords))
  cat(sprintf("  %s grouped/co-location leakage: %.1f%% of test points\n",
              if (f$group_leak_frac > 0) "[!]" else "[ok]", 100 * f$group_leak_frac))
  cat(sprintf("  %s CRS documented: %s\n",
              if (f$crs_documented) "[ok]" else "[!]", f$crs_documented))
  invisible(x)
}

#' Assemble a leakage report (scorecard)
#'
#' Combines a leakage diagnosis with optional optimism and validation
#' recommendation into a single printable scorecard for journal submission.
#'
#' @param diagnosis A [detect_leakage()] result.
#' @param optimism An optional [estimate_optimism()] result.
#' @param recommendation An optional [recommend_validation()] result.
#' @param deleak An optional [deleak_estimate()] result (the corrected accuracy).
#' @return An object of class `leakage_report`.
#' @export
report_leakage <- function(diagnosis, optimism = NULL, recommendation = NULL,
                           deleak = NULL) {
  stopifnot(inherits(diagnosis, "leakage_diagnosis"))
  structure(
    list(diagnosis = diagnosis, optimism = optimism, deleak = deleak,
         recommendation = recommendation, grade = .grade(diagnosis$SLI_rho)),
    class = "leakage_report")
}

#' @export
print.leakage_report <- function(x, ...) {
  cat("================ spLeakage report ================\n")
  cat(sprintf(" Leakage grade : %s\n", x$grade))
  cat("--------------------------------------------------\n")
  print(x$diagnosis)
  if (!is.null(x$optimism)) { cat("\n"); print(x$optimism) }
  if (!is.null(x$deleak)) { cat("\n"); print(x$deleak) }
  if (!is.null(x$recommendation)) { cat("\n"); print(x$recommendation) }
  cat("==================================================\n")
  invisible(x)
}
