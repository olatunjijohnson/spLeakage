make_st <- function(nloc = 40L, ntime = 10L, seed = 1L) {
  set.seed(seed)
  loc <- cbind(runif(nloc), runif(nloc)); g <- expand.grid(loc_id = seq_len(nloc), t = seq_len(ntime))
  data.frame(x = loc[g$loc_id, 1], y = loc[g$loc_id, 2], t = g$t, z = rnorm(nloc * ntime))
}

test_that("the joint space-time SLI catches leakage that space-only analysis misses", {
  d <- make_st()
  st <- detect_st_leakage(d, sample(rep_len(1:10, nrow(d))), time = "t", coords = c("x", "y"))
  expect_s3_class(st, "st_leakage")
  # all locations sampled -> space-only sees no leakage, but time/joint do.
  expect_lt(abs(st$sli_space), 0.05)
  expect_gt(st$sli_time, 0.1)
  expect_gt(st$sli_st, 0.05)
})

test_that("forward-chaining reduces the joint space-time leakage", {
  d <- make_st()
  st_rand <- detect_st_leakage(d, sample(rep_len(1:10, nrow(d))), time = "t", coords = c("x", "y"))
  st_fc   <- detect_st_leakage(d, temporal_kfold(d, k = 5, time = "t"), time = "t", coords = c("x", "y"))
  expect_lt(st_fc$sli_st, st_rand$sli_st)        # forward-chaining is less optimistic
})

test_that("supplying ranges and horizon is honoured", {
  d <- make_st()
  st <- detect_st_leakage(d, sample(rep_len(1:8, nrow(d))), time = "t", coords = c("x", "y"),
                          sp_range = 0.5, t_range = 3, horizon = 2)
  expect_equal(st$sp_range, 0.5); expect_equal(st$t_range, 3); expect_equal(st$horizon, 2)
})
