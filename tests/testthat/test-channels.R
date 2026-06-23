# --- feature-space leakage ------------------------------------------------------

test_that("feature-space leakage is positive when deployment extrapolates covariates", {
  set.seed(1)
  d <- data.frame(x = runif(200), y = runif(200), cov = runif(200, 0, 1))
  folds <- random_folds(200, k = 5)
  nd_extrap <- data.frame(cov = runif(300, 1.2, 2.2))   # outside sampled covariate range
  nd_same   <- data.frame(cov = runif(300, 0, 1))

  fl_x <- detect_feature_leakage(d, folds, covariates = "cov", newdata = nd_extrap)
  fl_s <- detect_feature_leakage(d, folds, covariates = "cov", newdata = nd_same)
  expect_s3_class(fl_x, "feature_leakage")
  expect_gt(fl_x$feature_sli, fl_s$feature_sli)   # extrapolation leaks more
  expect_gt(fl_x$feature_sli, 0)
})

test_that("feature leakage reports AOA without newdata", {
  set.seed(2)
  d <- data.frame(cov1 = rnorm(120), cov2 = rnorm(120))
  fl <- detect_feature_leakage(d, random_folds(120, 5), covariates = c("cov1", "cov2"))
  expect_false(fl$has_newdata)
  expect_true(fl$aoa_test >= 0 && fl$aoa_test <= 1)
})

# --- temporal leakage -----------------------------------------------------------

test_that("random split on time-stamped data has lookahead; temporal_kfold fixes it", {
  set.seed(3)
  d <- data.frame(x = runif(120), y = runif(120), t = 1:120)
  tl_rand <- detect_temporal_leakage(d, random_folds(120, 6), time = "t")
  expect_gt(tl_rand$lookahead_frac, 0.5)          # most test points trained on the future

  fc <- temporal_kfold(d, k = 6, time = "t")
  tl_fc <- detect_temporal_leakage(d, fc, time = "t")
  expect_equal(tl_fc$lookahead_frac, 0)           # forward-chaining: no lookahead
})

# --- pre-built fold lists flow through detect_leakage ---------------------------

test_that("detect_leakage accepts a pre-built fold list (e.g. temporal_kfold)", {
  d <- sim_spatial(120)
  d$t <- 1:120
  tgt <- prediction_target(grid = grid_target(), type = "grid")
  fc <- temporal_kfold(d, k = 5, time = "t")
  lk <- detect_leakage(d, fc, tgt, response = "z", coords = c("x", "y"))
  expect_s3_class(lk, "leakage_diagnosis")
})

# --- rsample adaptor ------------------------------------------------------------

test_that("rsample rset objects are accepted", {
  skip_if_not_installed("rsample")
  d <- sim_spatial(150)
  tgt <- prediction_target(grid = grid_target(), type = "grid")
  rs <- rsample::vfold_cv(d, v = 5)
  lk <- detect_leakage(d, rs, tgt, response = "z", coords = c("x", "y"))
  expect_s3_class(lk, "leakage_diagnosis")
  expect_equal(lk$n_folds, 5)
})
