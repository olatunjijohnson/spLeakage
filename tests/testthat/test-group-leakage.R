test_that("duplicated-location leakage is detected and fixed by group_kfold", {
  set.seed(1)
  base <- sim_spatial(n = 60L)
  d <- rbind(base, base[sample(60L, 15L), ])   # 75 rows, 60 unique locations
  folds <- sample(rep_len(1:5, nrow(d)))

  gl <- detect_group_leakage(d, folds, coords = c("x", "y"))
  expect_s3_class(gl, "group_leakage")
  expect_gt(gl$frac_leaked, 0)            # duplicates split across folds leak
  expect_gt(gl$n_split_groups, 0)
  expect_lt(gl$n_groups, nrow(d))         # fewer groups than rows (duplicates exist)

  gk <- group_kfold(d, k = 5L, coords = c("x", "y"))
  gl2 <- detect_group_leakage(d, gk, coords = c("x", "y"))
  expect_equal(gl2$frac_leaked, 0)        # group-aware folds eliminate it
})

test_that("no duplicates means no group leakage", {
  d <- sim_spatial(n = 80L)
  gl <- detect_group_leakage(d, random_folds(80L, k = 5L), coords = c("x", "y"))
  expect_equal(gl$frac_leaked, 0)
  expect_equal(gl$n_split_groups, 0)
  expect_equal(gl$n_multi_groups, 0)
})

test_that("explicit group column is honoured and fixable", {
  d <- sim_spatial(n = 60L); d$site <- rep(1:20, each = 3)   # 20 sites, 3 obs each
  folds <- random_folds(60L, k = 5L)
  gl <- detect_group_leakage(d, folds, group = "site", coords = c("x", "y"))
  expect_gt(gl$frac_leaked, 0)
  expect_equal(gl$n_groups, 20)

  gk <- group_kfold(d, k = 5L, group = "site")
  expect_equal(detect_group_leakage(d, gk, group = "site")$frac_leaked, 0)
})

test_that("near-duplicate grouping (tol > 0) merges jittered repeats", {
  set.seed(2)
  base <- sim_spatial(n = 40L)
  jit <- base[1:10, ]; jit$x <- jit$x + 1e-4; jit$y <- jit$y + 1e-4   # ~same place
  d <- rbind(base, jit)
  g_exact <- detect_group_leakage(d, random_folds(nrow(d), 5L), coords = c("x", "y"), tol = 0)
  g_near  <- detect_group_leakage(d, random_folds(nrow(d), 5L), coords = c("x", "y"), tol = 1e-3)
  expect_lt(g_near$n_groups, g_exact$n_groups)   # near-dup merges the jittered pairs
})

test_that("audit_workflow surfaces group leakage", {
  set.seed(3)
  base <- sim_spatial(n = 60L)
  d <- rbind(base, base[sample(60L, 12L), ])
  tgt <- prediction_target(grid = grid_target(), type = "grid")
  au <- audit_workflow(d, split = random_folds(nrow(d), 6L), target = tgt,
                       response = "z", coords = c("x", "y"))
  expect_true(au$flags$group_leak_frac > 0)
  expect_s3_class(au$group_leakage, "group_leakage")
})
