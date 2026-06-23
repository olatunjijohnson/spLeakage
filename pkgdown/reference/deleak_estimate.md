# De-leaked (bias-corrected) accuracy estimate

Reports the accuracy a model *would* achieve under a leakage-controlled,
design-matched validation – the de-leaked estimate practitioners should
report – with a fold (cluster) bootstrap confidence interval. Optionally
corrects a user-supplied `reported` metric value (e.g. a published RMSE)
using the estimated optimism ratio, without refitting the original
model.

## Usage

``` r
deleak_estimate(
  data,
  split,
  response,
  reported = NULL,
  target = NULL,
  predict_fun = NULL,
  metric = c("rmse", "mae", "brier", "logloss", "poisson"),
  control = c("block", "buffer"),
  coords = NULL,
  dependence = NULL,
  k_control = NULL,
  buffer = NULL,
  n_boot = 500L,
  level = 0.9
)
```

## Arguments

- data:

  An `sf` object, numeric coordinate matrix, or `data.frame`.

- split:

  The user's split: a `list(test=, train=)`, a list of test-index
  vectors, or a fold-id vector of length `nrow(data)`.

- response:

  Name of the numeric response column.

- reported:

  Optional reported metric value to correct (e.g. a published RMSE from
  a model not available for refitting). The de-leaked value is
  `reported * (controlled / user-CV)` error ratio.

- target:

  Optional
  [`prediction_target()`](https://olatunjijohnson.github.io/spLeakage/reference/prediction_target.md).
  When supplied, the controlled scheme is **target-matched** (the
  candidate scheme whose `SLI_rho` against the target is closest to
  zero) instead of a fixed spatial block – recommended when deployment
  is extrapolation, where a fixed block under-corrects.

- predict_fun:

  A function `(train_df, test_df)` returning numeric predictions for the
  test rows. Defaults to inverse-distance weighting on the coordinates.
  Custom learners (e.g. wrapping `lm`/`ranger`) let optimism be
  model-conditional.

- metric:

  Scoring rule (lower = better): `"rmse"` (default) or `"mae"` for
  continuous responses; `"brier"` or `"logloss"` for binary responses
  (predictions in `[0, 1]`); `"poisson"` deviance for counts. Use a
  proper rule for the response type so optimism is meaningful (e.g.
  Brier/log-loss for disease prevalence).

- control:

  The leakage-controlled comparison scheme: `"block"` (spatial k-means
  blocks, default) or `"buffer"` (buffer the user's folds by the
  practical range).

- coords:

  For non-`sf` input, the coordinate column names/indices.

- dependence:

  Optional
  [`estimate_dependence()`](https://olatunjijohnson.github.io/spLeakage/reference/estimate_dependence.md)
  object (used to set the buffer for `control = "buffer"`).

- k_control:

  Number of control folds (defaults to the number of user folds).

- buffer:

  Buffer distance for `control = "buffer"` (defaults to the practical
  range from `dependence`).

- n_boot:

  Number of fold-bootstrap resamples for the interval.

- level:

  Confidence level for the interval.

## Value

An object of class `deleak_estimate`.

## Examples

``` r
set.seed(1)
d <- data.frame(x = runif(150), y = runif(150), z = rnorm(150))
deleak_estimate(d, split = sample(rep_len(1:10, 150)), response = "z",
                coords = c("x", "y"), n_boot = 100)
#> <deleak_estimate>
#>   metric / control     : RMSE / block
#>   reported (your CV)   : 1.085
#>   de-leaked estimate   : 1.012   (90% CI 0.8921, 1.103)
#>   => reported accuracy inflated by -7% (error x0.93)
```
