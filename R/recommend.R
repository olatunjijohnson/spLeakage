# recommend_validation(): the design-aware recommendation engine (contribution C3).
# Estimand FIRST, then design x target. Design-basedness is ELICITED, never inferred
# from geometry; the clustering index is only a risk flag. See docs/VISION.md C3.

# Clark-Evans nearest-neighbour index (R < 1 clustered, ~1 random, > 1 regular).
# A risk flag only -- it detects spatial pattern, not whether the sample is a
# probability sample (which geometry cannot reveal).
.nn_index <- function(coords, info) {
  n <- nrow(coords)
  if (requireNamespace("FNN", quietly = TRUE) && !isTRUE(info$geographic)) {
    nn <- FNN::get.knn(coords, k = 1L)$nn.dist[, 1]
  } else {
    nn <- vapply(seq_len(n), function(i) {
      d <- sqrt(colSums((t(coords[-i, , drop = FALSE]) - coords[i, ])^2))
      min(d)
    }, numeric(1))
  }
  area <- prod(apply(coords, 2L, function(z) diff(range(z))))
  if (area <= 0) return(NA_real_)
  expected <- 0.5 / sqrt(n / area)
  mean(nn) / expected
}

# The decision matrix (VISION C3). Every branch is conditional and names what to
# avoid; design-based and model-based estimands are kept separate.
.decision <- function(estimand, design, target) {
  if (estimand == "population") {
    # Design-based map-accuracy estimand.
    if (design == "probability") {
      return(list(
        recommended = c("Design-based estimator (inclusion-weighted)", "Random CV"),
        avoid = "Spatial CV (introduces pessimistic bias for probability samples)",
        spatial = FALSE,
        rationale = paste("Population-mean map accuracy from a probability sample is",
                          "unbiased under design-based inference; spatial CV is not",
                          "appropriate here (Wadoux et al. 2021).")))
    }
    return(list(
      recommended = c("Design-based estimator with caution", "Model-assisted estimator"),
      avoid = "Naive random CV assuming it is unbiased",
      spatial = NA,
      rationale = paste("Population-mean accuracy from a non-probability sample lacks a",
                        "design basis; results are model-dependent. Treat with caution.")))
  }

  # Conditional predictive-skill estimand (the common ML case).
  if (design == "probability") {
    return(list(
      recommended = c("Random CV (unbiased for a probability sample)"),
      avoid = "Forcing spatial CV (over-pessimistic here)",
      spatial = FALSE,
      rationale = paste("For a probability sample, random CV gives unbiased predictive",
                        "skill; spatial CV would be pessimistic.")))
  }
  # clustered / convenience / unknown -> spatial CV matched to the target geometry.
  rec <- switch(target,
    interpolation = c("Buffered LOO CV", "NNDM LOO CV"),
    grid          = c("NNDM / kNNDM CV", "Spatial block CV"),
    newdata       = c("NNDM / kNNDM CV"))
  if (target == "newdata") rec <- c("Block CV + Area of Applicability", rec)
  list(
    recommended = rec,
    avoid = "Random k-fold CV (optimistic under spatial autocorrelation)",
    spatial = TRUE,
    rationale = paste0(
      "Conditional predictive skill from a ", design, " sample for a '", target,
      "' target: match the CV geometry to deployment (NNDM/kNNDM/buffered) so test",
      " points are as far from training as prediction points are from the sample."))
}

#' Recommend a validation strategy for a sampling design and prediction target
#'
#' Operationalises the design x estimand x target framework: it asks for the
#' estimand and the sampling design (which cannot be inferred from coordinates) and
#' returns a ranked, conditional recommendation -- including when spatial CV is
#' *not* appropriate. Any supplied data is used only to compute a clustering risk
#' flag, never to infer the design.
#'
#' @param data Optional `sf`/matrix/`data.frame`, used only for the clustering flag.
#' @param estimand `"prediction"` (conditional predictive skill at locations,
#'   default) or `"population"` (population-mean map accuracy over a region).
#' @param design `"unknown"` (default), `"probability"`, `"clustered"`, or
#'   `"convenience"`. This is an elicited fact about data collection.
#' @param target `"grid"` (wall-to-wall, default), `"interpolation"`, or `"newdata"`.
#' @param coords For non-`sf` input, the coordinate column names/indices.
#' @return An object of class `validation_recommendation`.
#' @examples
#' # A clustered sample to be mapped wall-to-wall: spatial CV is appropriate.
#' recommend_validation(estimand = "prediction", design = "clustered", target = "grid")
#' # A probability sample: random CV is correct, spatial CV is over-pessimistic.
#' recommend_validation(estimand = "prediction", design = "probability", target = "grid")
#' @export
recommend_validation <- function(data = NULL,
                                 estimand = c("prediction", "population"),
                                 design = c("unknown", "probability", "clustered", "convenience"),
                                 target = c("grid", "interpolation", "newdata"),
                                 coords = NULL) {
  estimand <- match.arg(estimand); design <- match.arg(design); target <- match.arg(target)
  flag <- NULL
  if (!is.null(data)) {
    xy <- .extract_coords(data, coords)
    R <- .nn_index(xy$coords, list(geographic = xy$geographic))
    flag <- list(nn_index = R, clustered = isTRUE(R < 0.85))
  }
  eff_design <- design
  if (design == "unknown") eff_design <- "clustered"  # conservative default
  rec <- .decision(estimand, eff_design, target)

  structure(
    list(estimand = estimand, design = design, target = target,
         recommended = rec$recommended, avoid = rec$avoid,
         spatial_cv_appropriate = rec$spatial, rationale = rec$rationale,
         clustering = flag, assumed_design = eff_design),
    class = "validation_recommendation")
}

#' @export
print.validation_recommendation <- function(x, ...) {
  cat("<validation_recommendation>\n")
  cat(sprintf("  estimand / design / target : %s / %s / %s\n",
              x$estimand, x$design, x$target))
  if (x$design == "unknown") {
    cat(sprintf("  NOTE: design unknown -> assuming '%s' (design cannot be inferred\n",
                x$assumed_design))
    cat("        from geometry; declare it for a definitive recommendation).\n")
  }
  sp <- if (is.na(x$spatial_cv_appropriate)) "depends (see rationale)"
        else if (x$spatial_cv_appropriate) "YES" else "NO"
  cat(sprintf("  spatial CV appropriate     : %s\n", sp))
  cat("  recommended:\n"); for (r in x$recommended) cat(sprintf("    - %s\n", r))
  cat(sprintf("  avoid: %s\n", x$avoid))
  if (!is.null(x$clustering) && is.finite(x$clustering$nn_index)) {
    cat(sprintf("  clustering flag (risk only): NN index = %.2f%s\n",
                x$clustering$nn_index,
                if (x$clustering$clustered) " (clustered)" else ""))
  }
  cat(sprintf("  rationale: %s\n", x$rationale))
  invisible(x)
}
