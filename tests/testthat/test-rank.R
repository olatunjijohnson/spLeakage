test_that("scheme ranking prefers a spatial scheme over random on clustered data", {
  d <- sim_spatial_clustered(n = 500L)
  tgt <- prediction_target(grid = grid_target(), type = "grid")
  rk <- rank_cv_schemes(d, tgt, response = "z", coords = c("x", "y"))
  expect_s3_class(rk, "scheme_ranking")
  expect_setequal(rk$ranking$scheme, c("random", "block", "buffered"))
  # random CV should be the most optimistic (largest SLI_rho), so NOT the best match.
  expect_false(identical(rk$best, "random"))
  rand_sli <- rk$ranking$sli_rho[rk$ranking$scheme == "random"]
  best_sli <- rk$ranking$abs_sli[1]
  expect_lt(best_sli, abs(rand_sli))
})

test_that("ranking is ordered by |SLI_rho| and names the best scheme", {
  d <- sim_spatial_clustered(n = 400L)
  tgt <- prediction_target(grid = grid_target(), type = "grid")
  rk <- rank_cv_schemes(d, tgt, response = "z", coords = c("x", "y"))
  expect_equal(rk$ranking$abs_sli, sort(rk$ranking$abs_sli))
  expect_equal(rk$best, rk$ranking$scheme[1])
})

test_that("a subset of schemes can be requested", {
  d <- sim_spatial(n = 200L)
  tgt <- prediction_target(grid = grid_target(), type = "grid")
  rk <- rank_cv_schemes(d, tgt, response = "z", coords = c("x", "y"),
                        schemes = c("random", "block"))
  expect_equal(nrow(rk$ranking), 2L)
})
