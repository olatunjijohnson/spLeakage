test_that("the shipped emulator loads and is well-formed", {
  em <- spLeakage:::.get_emulator()
  expect_s3_class(em, "optimism_emulator")
  expect_true(em$n_train > 200)
  expect_true(all(c("sli_rho", "range_rel", "signal", "nn_index") %in% em$features))
  expect_setequal(unique(em$model), c("idw", "rf", "gam"))
  expect_setequal(unique(em$response_type), c("gaussian", "poisson", "binomial"))
})

# A central, in-distribution feature vector for a given (model, response) cell:
# the column-wise median of that cell's (unstandardised) training features.
central_feats <- function(em, model = "idw", rt = "gaussian") {
  ix <- which(em$model == model & em$response_type == rt)
  X <- t(t(em$Xs[ix, , drop = FALSE]) * em$scale + em$center)
  stats::setNames(apply(X, 2L, stats::median), em$features)
}

test_that("an in-distribution feature vector gets an in-AOA finite estimate", {
  em <- spLeakage:::.get_emulator()
  res <- spLeakage:::.emu_predict_one(em, central_feats(em), "idw", "gaussian")
  expect_true(res$in_aoa)
  expect_true(is.finite(res$optimism_rel))
  expect_true(res$optimism_rel > -1 && res$optimism_rel < 1)
})

test_that("higher SLI predicts higher optimism (the C1 -> C2 relationship)", {
  em <- spLeakage:::.get_emulator()
  x <- central_feats(em)
  lo <- x; lo["sli_rho"] <- stats::quantile(
    (t(em$Xs) * em$scale + em$center)["sli_rho", ], 0.1)
  hi <- x; hi["sli_rho"] <- stats::quantile(
    (t(em$Xs) * em$scale + em$center)["sli_rho", ], 0.9)
  p_lo <- spLeakage:::.emu_predict_one(em, lo, "idw", "gaussian")$optimism_rel
  p_hi <- spLeakage:::.emu_predict_one(em, hi, "idw", "gaussian")$optimism_rel
  expect_gt(p_hi, p_lo)
})

test_that("emulator refuses an extreme out-of-AOA query", {
  em <- spLeakage:::.get_emulator()
  x <- central_feats(em); x["range_rel"] <- 50; x["sli_rho"] <- 0.9
  res <- spLeakage:::.emu_predict_one(em, x, "idw", "gaussian")
  expect_false(res$in_aoa)
})

test_that("predict_optimism runs end-to-end and returns a valid object", {
  d <- sim_gp(n = 300L, range = 0.2, signal = 0.6, design = "clustered")
  tgt <- prediction_target(grid = grid_target(), type = "grid")
  # May fall outside the AOA (different generator than the training study); the
  # AOA-refusal warning is expected behaviour, so allow it.
  pr <- suppressWarnings(predict_optimism(d, random_folds(nrow(d)), tgt,
                         response = "z", coords = c("x", "y"), model = "idw"))
  expect_s3_class(pr, "optimism_prediction")
  expect_equal(pr$response_type, "gaussian")
  expect_true(is.logical(pr$in_aoa))
  if (pr$in_aoa) expect_true(is.finite(pr$optimism_rel))
})
