test_that("estimate_dependence recovers a sensible range and signal", {
  d <- sim_spatial(n = 500L)
  dep <- estimate_dependence(d, response = "z", coords = c("x", "y"))
  expect_s3_class(dep, "sp_dependence")
  expect_gt(dep$practical_range, 0)
  expect_true(dep$signal_prop > 0.3)       # strong spatial signal in the sim
  expect_equal(dep$rho(0), 1)
  expect_lt(dep$rho(dep$practical_range), 0.1)  # ~0.05 at practical range
})

test_that("signed area equals the mean NN-distance gap (METHOD-SLI identity)", {
  d <- sim_spatial(n = 400L)
  tgt <- prediction_target(grid = grid_target(), type = "grid")
  lk <- detect_leakage(d, split = random_folds(nrow(d)), target = tgt,
                       response = "z", coords = c("x", "y"))
  gap <- mean(lk$fpred) - mean(lk$gobs[lk$tested])
  expect_equal(lk$A, gap, tolerance = 1e-8)
  expect_equal(lk$Aplus - lk$Aminus, lk$A, tolerance = 1e-8)
  expect_equal(lk$W, lk$Aplus + lk$Aminus, tolerance = 1e-8)
})

test_that("random split leaks more than a spatially blocked split", {
  d <- sim_spatial_clustered(n = 500L)
  xy <- as.matrix(d[, c("x", "y")])
  tgt <- prediction_target(grid = grid_target(), type = "grid")
  dep <- estimate_dependence(d, response = "z", coords = c("x", "y"))

  lk_rand <- detect_leakage(d, split = random_folds(nrow(d)), target = tgt,
                            dependence = dep, coords = c("x", "y"))
  lk_block <- detect_leakage(d, split = block_folds(xy), target = tgt,
                             dependence = dep, coords = c("x", "y"))

  # Random CV is optimistic; blocking reduces (or reverses) the leakage.
  expect_gt(lk_rand$SLI_rho, lk_block$SLI_rho)
  expect_gt(lk_rand$SLI_rho, 0)
  expect_gt(lk_rand$SLI_d, lk_block$SLI_d)
})

test_that("sli() accessor and signs behave", {
  d <- sim_spatial(n = 300L)
  tgt <- prediction_target(grid = grid_target(), type = "grid")
  lk <- detect_leakage(d, split = random_folds(nrow(d)), target = tgt,
                       response = "z", coords = c("x", "y"))
  expect_equal(sli(lk), lk$SLI_rho)
  expect_equal(sli(lk, "d"), lk$SLI_d)
  expect_true(lk$SLI_rho >= -1 && lk$SLI_rho <= 1)
  expect_equal(unname(mean(lk$leak_point[lk$tested])), lk$SLI_rho, tolerance = 1e-10)
})

test_that("split parser accepts folds list, fold-ids, and single split", {
  d <- sim_spatial(n = 120L)
  tgt <- prediction_target(grid = grid_target(10L), type = "grid")
  dep <- estimate_dependence(d, response = "z", coords = c("x", "y"))
  ids <- random_folds(nrow(d), k = 5L)
  fold_list <- split(seq_len(nrow(d)), ids)

  l1 <- detect_leakage(d, split = ids, target = tgt, dependence = dep, coords = c("x", "y"))
  l2 <- detect_leakage(d, split = fold_list, target = tgt, dependence = dep, coords = c("x", "y"))
  expect_equal(l1$SLI_rho, l2$SLI_rho, tolerance = 1e-10)

  test_idx <- which(ids == 1L)
  l3 <- detect_leakage(d, split = list(test = test_idx), target = tgt,
                       dependence = dep, coords = c("x", "y"))
  expect_equal(l3$n_test, length(test_idx))
})

test_that("works with sf input and geographic CRS", {
  skip_if_not_installed("sf")
  d <- sim_spatial(n = 200L)
  # Map unit square into a small lon/lat patch.
  d$lon <- -2 + d$x; d$lat <- 53 + d$y
  sfd <- sf::st_as_sf(d, coords = c("lon", "lat"), crs = 4326)
  gt <- grid_target(12L)
  gsf <- sf::st_as_sf(data.frame(lon = -2 + gt[, 1], lat = 53 + gt[, 2]),
                      coords = c("lon", "lat"), crs = 4326)
  tgt <- prediction_target(grid = gsf, type = "grid")
  lk <- detect_leakage(sfd, split = random_folds(nrow(d)), target = tgt,
                       response = "z")
  expect_s3_class(lk, "leakage_diagnosis")
  expect_true(is.finite(lk$SLI_rho))
})
