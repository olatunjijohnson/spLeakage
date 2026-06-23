test_that("trend_strength is high for a trend-dominated field, low for pure noise", {
  set.seed(1)
  trendy <- data.frame(x = runif(200), y = runif(200))
  trendy$z <- 3 * trendy$x - 2 * trendy$y + rnorm(200, 0, 0.2)   # strong linear trend
  noise <- data.frame(x = runif(200), y = runif(200), z = rnorm(200))

  st <- trend_strength(trendy, "z", coords = c("x", "y"))
  sn <- trend_strength(noise, "z", coords = c("x", "y"))
  expect_true(st > 0.8)
  expect_true(sn < 0.2)
  expect_gt(st, sn)
})

test_that("trend_strength returns a proportion in [0,1]", {
  d <- sim_spatial(n = 150L)
  s <- trend_strength(d, "z", coords = c("x", "y"))
  expect_true(s >= 0 && s <= 1)
})

test_that("audit_workflow reports the trend channel", {
  set.seed(2)
  d <- data.frame(x = runif(200), y = runif(200)); d$z <- 4 * d$x + rnorm(200, 0, 0.3)
  tgt <- prediction_target(grid = grid_target(), type = "grid")
  au <- audit_workflow(d, split = random_folds(200, 6), target = tgt,
                       response = "z", coords = c("x", "y"))
  expect_true(is.finite(au$flags$trend_strength))
  expect_gt(au$flags$trend_strength, 0.5)            # this field is trend-dominated
  expect_output(print(au), "trend strength")
})
