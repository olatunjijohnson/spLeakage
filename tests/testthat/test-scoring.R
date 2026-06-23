# Binary (prevalence-like) and count fields for proper scoring rules.
sim_binary <- function(n = 400L, seed = 1L) {
  set.seed(seed)
  nc <- 8; ce <- cbind(runif(nc, .1, .9), runif(nc, .1, .9)); cl <- sample(nc, n, TRUE)
  xy <- ce[cl, ] + cbind(rnorm(n, 0, .04), rnorm(n, 0, .04))
  xy <- cbind(pmin(pmax(xy[, 1], 0), 1), pmin(pmax(xy[, 2], 0), 1))
  p <- plogis(2 * sin(2 * pi * xy[, 1]) + 2 * cos(2 * pi * xy[, 2]))
  data.frame(x = xy[, 1], y = xy[, 2], z = rbinom(n, 1, p))
}

test_that("proper scoring rules run for binary responses (Brier, log-loss)", {
  d <- sim_binary()
  for (mt in c("brier", "logloss")) {
    opt <- estimate_optimism(d, random_folds(nrow(d)), response = "z",
                             coords = c("x", "y"), metric = mt)
    expect_true(is.finite(opt$E_cv) && is.finite(opt$E_control))
    expect_gt(opt$optimism_rel, 0)             # random CV optimistic on clustered binary data
  }
})

test_that("Poisson deviance runs for counts", {
  set.seed(2)
  d <- sim_binary(); d$z <- rpois(nrow(d), exp(0.5 + d$x))
  opt <- estimate_optimism(d, random_folds(nrow(d)), response = "z",
                           coords = c("x", "y"), metric = "poisson")
  expect_true(is.finite(opt$E_cv))
})

test_that("deleak_estimate supports proper scoring rules", {
  d <- sim_binary()
  dl <- deleak_estimate(d, random_folds(nrow(d)), response = "z", coords = c("x", "y"),
                        metric = "brier", n_boot = 100L)
  expect_equal(dl$metric, "brier")
  expect_gt(dl$deleaked, dl$reported_cv)
  expect_length(dl$deleaked_ci, 2L)
})

test_that("the metric specs are proper / sane", {
  # perfect prediction scores zero; worse predictions score higher.
  brier <- spLeakage:::.metric_fun("brier")
  ll <- spLeakage:::.metric_fun("logloss")
  expect_equal(brier(c(1, 0), c(1, 0)), 0)
  expect_lt(brier(c(1, 0), c(0.9, 0.1)), brier(c(1, 0), c(0.6, 0.4)))
  expect_lt(ll(c(1, 0), c(0.9, 0.1)), ll(c(1, 0), c(0.6, 0.4)))
})
