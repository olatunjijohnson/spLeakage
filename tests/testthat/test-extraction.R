test_that("extraction-overlap leakage is detected and removed by a 2r buffer", {
  set.seed(1)
  d <- data.frame(x = runif(300), y = runif(300), z = rnorm(300))
  folds <- random_folds(300, k = 5)
  r <- 0.08

  el <- detect_extraction_leakage(d, folds, radius = r, coords = c("x", "y"))
  expect_s3_class(el, "extraction_leakage")
  expect_gt(el$frac_overlap, 0)                 # dense random data -> windows overlap
  expect_gte(el$frac_overlap, el$frac_contains) # overlap (<2r) >= contains (<r)
  expect_equal(el$recommended_buffer, 2 * r)

  # Buffer the folds by 2r -> no remaining extraction overlap.
  xy <- as.matrix(d[, c("x", "y")])
  base <- lapply(split(seq_len(300), folds), function(i) list(test = i, train = setdiff(seq_len(300), i)))
  buffered <- spLeakage:::.buffer_folds(xy, base, 2 * r, list(geographic = FALSE))
  el2 <- detect_extraction_leakage(d, buffered, radius = r, coords = c("x", "y"))
  expect_equal(el2$frac_overlap, 0)
})

test_that("a tiny extraction radius yields little overlap", {
  set.seed(2)
  d <- data.frame(x = runif(300), y = runif(300), z = rnorm(300))
  el <- detect_extraction_leakage(d, random_folds(300, 5), radius = 0.005, coords = c("x", "y"))
  expect_lt(el$frac_overlap, 0.2)
})

test_that("radius must be positive", {
  d <- data.frame(x = runif(50), y = runif(50), z = rnorm(50))
  expect_error(detect_extraction_leakage(d, random_folds(50, 5), radius = -1, coords = c("x", "y")),
               "positive")
})
