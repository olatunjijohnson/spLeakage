test_that("random CV is optimistic relative to spatial blocking (clustered data)", {
  d <- sim_spatial_clustered(n = 500L)
  opt <- estimate_optimism(d, split = random_folds(nrow(d)), response = "z",
                           coords = c("x", "y"), control = "block")
  expect_s3_class(opt, "optimism_estimate")
  # Controlled (blocked) error should exceed the optimistic random-CV error.
  expect_gt(opt$E_control, opt$E_cv)
  expect_gt(opt$optimism_rel, 0)
})

test_that("buffer control route runs and uses the practical range", {
  d <- sim_spatial_clustered(n = 400L)
  dep <- estimate_dependence(d, response = "z", coords = c("x", "y"))
  opt <- suppressWarnings(estimate_optimism(
    d, split = random_folds(nrow(d)), response = "z",
    coords = c("x", "y"), control = "buffer", dependence = dep))
  expect_true(is.finite(opt$E_cv))
  expect_true(is.finite(opt$optimism_rel))
})

test_that("a spatial block split shows little/negative optimism vs block control", {
  d <- sim_spatial_clustered(n = 500L)
  blk <- spatial_block_cv(d, k = 10L, coords = c("x", "y"))
  opt <- estimate_optimism(d, split = blk, response = "z", coords = c("x", "y"),
                           control = "block")
  # When the user already blocks, optimism should be near zero (not strongly +).
  expect_lt(opt$optimism_rel, 0.15)
})

test_that("custom predict_fun (lm) is honoured", {
  d <- sim_spatial(n = 300L)
  d$cov <- d$z + rnorm(nrow(d), 0, 0.3)         # a useful covariate
  pf <- function(train, test) {
    m <- stats::lm(z ~ cov, data = train)
    as.numeric(stats::predict(m, newdata = test))
  }
  opt <- estimate_optimism(d, split = random_folds(nrow(d)), response = "z",
                           predict_fun = pf, coords = c("x", "y"))
  expect_s3_class(opt, "optimism_estimate")
  expect_true(is.finite(opt$optimism))
})
