test_that("design_validation returns a usable, budgeted, weighted design", {
  set.seed(1)
  d <- sim_spatial_clustered(n = 300L)
  tgt <- prediction_target(grid = grid_target(20L), type = "grid")
  des <- design_validation(d, tgt, budget = 30L, response = "z", coords = c("x", "y"),
                           allocation = "optimal")
  expect_s3_class(des, "validation_design")
  expect_true(abs(des$n_selected - 30L) <= 6L)         # ~budget (rounding across strata)
  expect_equal(nrow(des$locations), des$n_selected)
  expect_equal(sum(des$weights), 1, tolerance = 1e-8)
  expect_gt(des$se_reduction_vs_random, 0)             # optimal beats random
})

test_that("optimal allocation puts more points where the error is larger", {
  set.seed(2)
  d <- sim_spatial_clustered(n = 300L)
  tgt <- prediction_target(grid = grid_target(20L), type = "grid")
  des <- design_validation(d, tgt, budget = 40L, response = "z", coords = c("x", "y"),
                           allocation = "optimal")
  s <- des$strata
  # the highest-error stratum gets more points than the lowest-error stratum.
  expect_gte(s$n[which.max(s$mean_error)], s$n[which.min(s$mean_error)])
})

test_that("random allocation reports no SE gain", {
  set.seed(3)
  d <- sim_spatial(n = 200L)
  tgt <- prediction_target(grid = grid_target(), type = "grid")
  des <- design_validation(d, tgt, budget = 25L, response = "z", coords = c("x", "y"),
                           allocation = "random")
  expect_equal(des$allocation, "random")
  expect_equal(des$n_selected, 25L)
})
