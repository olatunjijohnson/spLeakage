# Grouped / duplicated-location leakage (contribution C4). A discrete, exact form of
# leakage complementing the distance-based SLI: when repeated observations at the
# same site (or members of the same group/cluster) are split across train and test,
# the test set is not independent regardless of the autocorrelation range. Surfaced
# as the dominant real-world channel in docs/CASE-STUDY-NIGERIA.md.

# Group id per observation, from an explicit grouping or from coordinates (exact
# duplicates when tol = 0; connected components within `tol` distance otherwise).
.group_from_coords <- function(coords, info, tol = 0) {
  n <- nrow(coords)
  if (tol <= 0) {
    key <- paste(round(coords[, 1], 8), round(coords[, 2], 8), sep = "_")
    return(match(key, unique(key)))
  }
  parent <- seq_len(n)
  find <- function(x) { while (parent[x] != x) x <- parent[x]; x }   # read-only
  for (i in seq_len(n - 1L)) {
    d <- .nn_dist(coords[(i + 1L):n, , drop = FALSE], coords[i, , drop = FALSE], info)
    for (j in which(d <= tol)) {
      a <- find(i); b <- find(i + j)
      if (a != b) parent[b] <- a       # union in the local scope
    }
  }
  roots <- vapply(seq_len(n), find, integer(1))
  match(roots, unique(roots))
}

.resolve_groups <- function(data, group, coords, tol) {
  if (!is.null(group)) {
    df <- if (inherits(data, "sf")) sf::st_drop_geometry(data) else as.data.frame(data)
    if (!group %in% names(df)) stop(sprintf("Group column '%s' not found.", group))
    return(list(grp = as.integer(factor(df[[group]])), source = group))
  }
  xy <- .extract_coords(data, coords)
  list(grp = .group_from_coords(xy$coords, list(crs = xy$crs, geographic = xy$geographic), tol),
       source = sprintf("coordinates (tol = %g)", tol))
}

#' Detect grouped / duplicated-location leakage in a split
#'
#' Flags test observations that share a group (e.g. the same site, repeated survey,
#' household, or plot) with a training observation in the same fold -- exact leakage
#' that the distance-based [detect_leakage()] index registers only at distance zero.
#'
#' @param data An `sf` object, numeric matrix, or `data.frame`.
#' @param split A split specification (see [detect_leakage()]).
#' @param group Optional name of a grouping column. If omitted, groups are derived
#'   from the coordinates.
#' @param coords For non-`sf` input, the coordinate column names/indices.
#' @param tol When grouping by coordinates, the distance within which points are
#'   treated as the same location (`0` = exact duplicates).
#' @return An object of class `group_leakage`.
#' @examples
#' # Two sites, each measured twice; a random split co-locates train and test.
#' d <- data.frame(x = c(1, 1, 2, 2), y = c(1, 1, 2, 2), z = rnorm(4))
#' detect_group_leakage(d, split = c(1, 2, 1, 2), coords = c("x", "y"))
#' @export
detect_group_leakage <- function(data, split, group = NULL, coords = NULL, tol = 0) {
  rg <- .resolve_groups(data, group, coords, tol)
  grp <- rg$grp; n <- length(grp)
  folds <- .parse_split(split, n)
  leaked <- logical(n); fold_id <- rep(NA_integer_, n); split_groups <- integer(0)
  for (k in seq_along(folds)) {
    f <- folds[[k]]
    if (!length(f$test) || !length(f$train)) next
    tr_groups <- unique(grp[f$train])
    hit <- grp[f$test] %in% tr_groups
    leaked[f$test] <- hit; fold_id[f$test] <- k
    split_groups <- c(split_groups, unique(grp[f$test][hit]))
  }
  tested <- which(!is.na(fold_id))
  gs <- table(grp)
  structure(
    list(n = n, n_test = length(tested), n_leaked = sum(leaked[tested]),
         frac_leaked = if (length(tested)) mean(leaked[tested]) else 0,
         n_groups = length(gs), n_multi_groups = sum(gs > 1L),
         n_split_groups = length(unique(split_groups)),
         group = grp, leaked = leaked, fold_id = fold_id,
         tol = tol, source = rg$source),
    class = "group_leakage")
}

#' Group-aware k-fold assignment (the fix for grouped leakage)
#'
#' Assigns whole groups to folds so that members of a group (e.g. co-located
#' records) are never split across train and test. The remedy for the leakage that
#' [detect_group_leakage()] flags.
#'
#' @inheritParams detect_group_leakage
#' @param k Number of folds.
#' @return An integer fold-id vector of length `nrow(data)`.
#' @export
group_kfold <- function(data, k = 10L, group = NULL, coords = NULL, tol = 0) {
  grp <- .resolve_groups(data, group, coords, tol)$grp
  sizes <- table(grp)
  k <- min(k, length(sizes))
  ord <- names(sizes)[order(-as.integer(sizes))]   # largest groups first
  load <- numeric(k); gfold <- stats::setNames(integer(length(sizes)), names(sizes))
  for (g in ord) { f <- which.min(load); gfold[g] <- f; load[f] <- load[f] + sizes[g] }
  as.integer(gfold[as.character(grp)])
}

#' @export
print.group_leakage <- function(x, ...) {
  cat("<group_leakage>\n")
  cat(sprintf("  grouping        : %s  |  n = %d, groups = %d, multi-member = %d\n",
              x$source, x$n, x$n_groups, x$n_multi_groups))
  cat(sprintf("  test leaked via shared group : %d / %d (%.1f%%)\n",
              x$n_leaked, x$n_test, 100 * x$frac_leaked))
  cat(sprintf("  groups split across folds    : %d\n", x$n_split_groups))
  if (x$frac_leaked > 0) {
    cat("  [!] fix: group_kfold() keeps co-located/grouped records together\n")
  }
  invisible(x)
}
