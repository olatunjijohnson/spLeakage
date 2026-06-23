# S3 methods and accessors for leakage_diagnosis objects.

#' Extract the Spatial Leakage Index
#'
#' @param x A `leakage_diagnosis` object.
#' @param type `"rho"` (dependence form, default) or `"d"` (distance form).
#' @return A single numeric SLI value.
#' @export
sli <- function(x, type = c("rho", "d")) {
  stopifnot(inherits(x, "leakage_diagnosis"))
  type <- match.arg(type)
  if (type == "rho") x$SLI_rho else x$SLI_d
}

# Human-readable verdict from the signed index.
.leak_verdict <- function(sli_rho, delta) {
  dir <- if (sli_rho > 0.02) "OPTIMISTIC leakage"
         else if (sli_rho < -0.02) "PESSIMISTIC (anti-leakage)"
         else "well matched"
  cross <- if (abs(delta) < 0.5) " (ECDFs cross; inspect plot)" else ""
  paste0(dir, cross)
}

#' @export
print.leakage_diagnosis <- function(x, ...) {
  cat("<leakage_diagnosis>\n")
  cat(sprintf("  target            : %s   |  n = %d, test = %d, folds = %d%s\n",
              x$target_type, x$n, x$n_test, x$n_folds,
              if (isTRUE(x$anisotropic)) "  [anisotropic]" else ""))
  ci_r <- if (!is.null(x$SLI_rho_ci))
    sprintf("  90%% CI [%+.3f, %+.3f]", x$SLI_rho_ci[1], x$SLI_rho_ci[2]) else ""
  ci_d <- if (!is.null(x$SLI_d_ci))
    sprintf("  90%% CI [%+.3f, %+.3f]", x$SLI_d_ci[1], x$SLI_d_ci[2]) else ""
  cat(sprintf("  SLI_rho (signed)  : %+.3f   [%s]%s\n",
              x$SLI_rho, .leak_verdict(x$SLI_rho, x$delta), ci_r))
  cat(sprintf("  SLI_d  (signed)   : %+.3f   (A = %+.4g, phi = %.4g)%s\n",
              x$SLI_d, x$A, x$phi, ci_d))
  cat(sprintf("  retained corr.    : c_obs = %.3f vs c_pred = %.3f\n",
              x$c_obs, x$c_pred))
  cat(sprintf("  W (NNDM) / delta  : %.4g / %+.2f\n", x$W, x$delta))
  invisible(x)
}

#' @export
summary.leakage_diagnosis <- function(object, ...) {
  print(object)
  if (length(object$fold_sli) > 1L) {
    cat("\n  per-fold leakage (mean rho excess):\n")
    fs <- object$fold_sli
    for (k in seq_along(fs)) {
      cat(sprintf("    fold %-3s : %+.3f\n", names(fs)[k], fs[k]))
    }
  }
  invisible(object)
}

#' Plot a leakage diagnosis
#'
#' @param x A `leakage_diagnosis` object.
#' @param which `"ecdf"` (observed vs target NN-distance ECDFs) or `"map"`
#'   (per-point leakage contribution in space).
#' @param ... Passed to the underlying plot call.
#' @return `x`, invisibly.
#' @export
plot.leakage_diagnosis <- function(x, which = c("ecdf", "map"), ...) {
  which <- match.arg(which)
  if (which == "ecdf") {
    e <- x$ecdf
    plot(e$r, e$Gobs, type = "s", col = "#d1495b", lwd = 2,
         xlab = "nearest-neighbour distance", ylab = "ECDF",
         main = "Observed (test->train) vs target (pred->sample)",
         ylim = c(0, 1), ...)
    graphics::lines(e$r, e$Gpred, type = "s", col = "#00798c", lwd = 2)
    graphics::abline(v = x$phi, lty = 3, col = "grey50")
    graphics::legend("bottomright", bty = "n",
                     col = c("#d1495b", "#00798c", "grey50"),
                     lwd = c(2, 2, 1), lty = c(1, 1, 3),
                     legend = c("observed (CV)", "target (deployment)", "phi"))
  } else {
    co <- x$coords; lp <- x$leak_point
    keep <- !is.na(lp)
    rng <- max(abs(lp[keep]))
    pal <- grDevices::colorRampPalette(c("#00798c", "grey90", "#d1495b"))(100)
    idx <- as.integer(cut(lp[keep], breaks = seq(-rng, rng, length.out = 101),
                          include.lowest = TRUE))
    plot(co[keep, 1], co[keep, 2], col = pal[idx], pch = 19,
         xlab = "x", ylab = "y", main = "Per-point leakage contribution", ...)
  }
  invisible(x)
}
