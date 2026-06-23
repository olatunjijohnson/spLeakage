# Anisotropic field: correlation reaches much farther along the x-axis than y.
sim_anisotropic <- function(n = 500L, seed = 7L) {
  set.seed(seed)
  xy <- cbind(runif(n), runif(n))
  # Low frequency in x (long range), high frequency in y (short range).
  field <- sin(2 * pi * xy[, 1]) + sin(8 * pi * xy[, 2])
  data.frame(x = xy[, 1], y = xy[, 2], z = field + rnorm(n, 0, 0.1))
}

test_that("isotropy reduces to the rotation-invariant baseline", {
  d <- sim_spatial(n = 300L)
  tgt <- prediction_target(grid = grid_target(), type = "grid")
  dep_iso <- estimate_dependence(d, response = "z", coords = c("x", "y"))
  # ratio = 1 is a pure rotation; distances (hence SLI) must be unchanged.
  dep_rot <- estimate_dependence(d, response = "z", coords = c("x", "y"),
                                 anisotropy = list(angle = pi / 3, ratio = 1))
  l_iso <- detect_leakage(d, random_folds(nrow(d)), tgt, dependence = dep_iso, coords = c("x", "y"))
  l_rot <- detect_leakage(d, random_folds(nrow(d)), tgt, dependence = dep_rot, coords = c("x", "y"))
  expect_equal(l_iso$SLI_rho, l_rot$SLI_rho, tolerance = 1e-8)
})

test_that("auto anisotropy is detected on an anisotropic field", {
  d <- sim_anisotropic(n = 600L)
  dep <- estimate_dependence(d, response = "z", coords = c("x", "y"),
                             anisotropy = "auto", n_bins = 12L)
  expect_false(is.null(dep$anisotropy))
  # Minor/major ratio should be clearly < 1 (strong directionality).
  expect_lt(dep$anisotropy$ratio, 0.9)
})

test_that("anisotropy-aware SLI differs from the isotropic SLI", {
  d <- sim_anisotropic(n = 600L)
  tgt <- prediction_target(grid = grid_target(), type = "grid")
  dep_iso <- estimate_dependence(d, response = "z", coords = c("x", "y"), n_bins = 12L)
  dep_ani <- estimate_dependence(d, response = "z", coords = c("x", "y"),
                                 anisotropy = "auto", n_bins = 12L)
  l_iso <- detect_leakage(d, random_folds(nrow(d)), tgt, dependence = dep_iso, coords = c("x", "y"))
  l_ani <- detect_leakage(d, random_folds(nrow(d)), tgt, dependence = dep_ani, coords = c("x", "y"))
  expect_true(l_ani$anisotropic)
  expect_false(isTRUE(all.equal(l_iso$SLI_rho, l_ani$SLI_rho)))
})

test_that("Monte-Carlo uncertainty brackets the point estimate", {
  d <- sim_spatial_clustered(n = 500L)
  tgt <- prediction_target(grid = grid_target(), type = "grid")
  lk <- detect_leakage(d, random_folds(nrow(d)), tgt, response = "z",
                       coords = c("x", "y"), n_boot = 200L)
  expect_length(lk$SLI_rho_ci, 2L)
  expect_lt(lk$SLI_rho_ci[1], lk$SLI_rho_ci[2])
  expect_gte(lk$SLI_rho, lk$SLI_rho_ci[1] - 1e-6)
  expect_lte(lk$SLI_rho, lk$SLI_rho_ci[2] + 1e-6)
  expect_length(lk$boot$SLI_rho, 200L)
})

test_that("geographic anisotropy request warns and is ignored", {
  d <- sim_spatial(n = 150L)
  d$lon <- -2 + d$x; d$lat <- 53 + d$y
  sfd <- sf::st_as_sf(d, coords = c("lon", "lat"), crs = 4326)
  expect_warning(
    dep <- estimate_dependence(sfd, response = "z", anisotropy = "auto"),
    "not supported for geographic")
  expect_null(dep$anisotropy)
})
