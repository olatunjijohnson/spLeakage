test_that("probability sample is steered AWAY from spatial CV", {
  rec <- recommend_validation(estimand = "prediction", design = "probability",
                              target = "grid")
  expect_false(rec$spatial_cv_appropriate)
  expect_match(rec$avoid, "spatial", ignore.case = TRUE)
})

test_that("clustered sample for a map is steered TOWARD spatial CV", {
  rec <- recommend_validation(estimand = "prediction", design = "clustered",
                              target = "grid")
  expect_true(rec$spatial_cv_appropriate)
  expect_match(rec$avoid, "Random", ignore.case = TRUE)
})

test_that("population estimand + probability sample -> design-based, not spatial", {
  rec <- recommend_validation(estimand = "population", design = "probability",
                              target = "grid")
  expect_false(rec$spatial_cv_appropriate)
  expect_match(paste(rec$recommended, collapse = " "), "[Dd]esign-based")
})

test_that("unknown design is conservative and flagged, not silently inferred", {
  rec <- recommend_validation(estimand = "prediction", design = "unknown",
                              target = "grid")
  expect_equal(rec$assumed_design, "clustered")
})

test_that("clustering flag is computed from data but does not set the design", {
  d <- sim_spatial_clustered(n = 400L)
  rec <- recommend_validation(d, estimand = "prediction", design = "unknown",
                              target = "grid", coords = c("x", "y"))
  expect_true(rec$clustering$clustered)       # geometry IS clustered
  expect_equal(rec$design, "unknown")         # but design stays unknown
})

test_that("audit_workflow and report_leakage assemble a scorecard", {
  d <- sim_spatial_clustered(n = 400L)
  tgt <- prediction_target(grid = grid_target(), type = "grid")
  au <- audit_workflow(d, split = random_folds(nrow(d)), target = tgt,
                       response = "z", coords = c("x", "y"))
  expect_s3_class(au, "workflow_audit")
  expect_true(au$flags$leakage_grade %in% c("A", "B", "C", "D", "F"))

  lk <- detect_leakage(d, random_folds(nrow(d)), tgt, response = "z", coords = c("x", "y"))
  opt <- estimate_optimism(d, random_folds(nrow(d)), response = "z", coords = c("x", "y"))
  rec <- recommend_validation(estimand = "prediction", design = "clustered", target = "grid")
  rep <- report_leakage(lk, optimism = opt, recommendation = rec)
  expect_s3_class(rep, "leakage_report")
  expect_output(print(rep), "spLeakage report")
})
