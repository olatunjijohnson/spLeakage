test_that("de-leaked estimate exceeds the reported CV error on clustered data", {
  d <- sim_spatial_clustered(n = 500L)
  dl <- deleak_estimate(d, split = random_folds(nrow(d)), response = "z",
                        coords = c("x", "y"), n_boot = 200L)
  expect_s3_class(dl, "deleak_estimate")
  expect_gt(dl$deleaked, dl$reported_cv)        # correction inflates the error
  expect_gt(dl$optimism_rel, 0)
  expect_gt(dl$ratio, 1)
})

test_that("the de-leaked confidence interval brackets the point estimate", {
  d <- sim_spatial_clustered(n = 400L)
  dl <- deleak_estimate(d, split = random_folds(nrow(d)), response = "z",
                        coords = c("x", "y"), n_boot = 300L)
  expect_length(dl$deleaked_ci, 2L)
  expect_lt(dl$deleaked_ci[1], dl$deleaked_ci[2])
  expect_gte(dl$deleaked, dl$deleaked_ci[1] - 1e-6)
  expect_lte(dl$deleaked, dl$deleaked_ci[2] + 1e-6)
})

test_that("a reported value is corrected by the optimism ratio (no refit needed)", {
  d <- sim_spatial_clustered(n = 400L)
  dl <- deleak_estimate(d, split = random_folds(nrow(d)), response = "z",
                        coords = c("x", "y"), reported = 0.5, n_boot = 100L)
  expect_equal(dl$anchored, 0.5 * dl$ratio, tolerance = 1e-8)
  expect_gt(dl$anchored, 0.5)                    # published value was optimistic
  expect_length(dl$anchored_ci, 2L)
})

test_that("target-aware de-leak matches the declared target better than fixed block", {
  d <- sim_spatial_clustered(n = 500L)
  tgt <- prediction_target(grid = grid_target(), type = "grid")
  dl_t <- deleak_estimate(d, split = random_folds(nrow(d)), response = "z", target = tgt,
                          coords = c("x", "y"), n_boot = 100L)
  expect_match(dl_t$control, "target-matched")
  expect_true(is.finite(dl_t$matched_sli))
  # The matched control's residual leakage is small (near zero) by construction.
  expect_lt(abs(dl_t$matched_sli), 0.1)
  expect_gt(dl_t$deleaked, dl_t$reported_cv)
})

test_that("a spatially blocked split needs little correction", {
  d <- sim_spatial_clustered(n = 500L)
  blk <- spatial_block_cv(d, k = 10L, coords = c("x", "y"))
  dl <- deleak_estimate(d, split = blk, response = "z", coords = c("x", "y"),
                        n_boot = 150L)
  expect_lt(dl$optimism_rel, 0.15)              # already de-leaked
})
