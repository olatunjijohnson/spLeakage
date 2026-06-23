test_that("k = 1 reproduces the default single-neighbour SLI", {
  d <- sim_spatial_clustered(n = 400L)
  tgt <- prediction_target(grid = grid_target(), type = "grid")
  dep <- estimate_dependence(d, response = "z", coords = c("x", "y"))
  l_def <- detect_leakage(d, random_folds(nrow(d)), tgt, dependence = dep, coords = c("x", "y"))
  l_k1  <- detect_leakage(d, random_folds(nrow(d)), tgt, dependence = dep, coords = c("x", "y"), k = 1L)
  expect_equal(l_def$SLI_rho, l_k1$SLI_rho, tolerance = 1e-10)
})

test_that("k > 1 runs and changes the dependence-form SLI", {
  d <- sim_spatial_clustered(n = 400L)
  tgt <- prediction_target(grid = grid_target(), type = "grid")
  dep <- estimate_dependence(d, response = "z", coords = c("x", "y"))
  l1 <- detect_leakage(d, random_folds(nrow(d)), tgt, dependence = dep, coords = c("x", "y"), k = 1L)
  l5 <- detect_leakage(d, random_folds(nrow(d)), tgt, dependence = dep, coords = c("x", "y"), k = 5L)
  expect_true(is.finite(l5$SLI_rho))
  expect_false(isTRUE(all.equal(l1$SLI_rho, l5$SLI_rho)))   # k>1 uses more neighbours
  expect_equal(l1$SLI_d, l5$SLI_d, tolerance = 1e-10)        # distance form unchanged by k
})
