# Estimate the optimism (accuracy inflation) of a cross-validation scheme

Re-evaluates a model under the user's split and under a
leakage-controlled scheme, and reports the gap as optimism: how much the
user's reported accuracy is inflated. Positive optimism means the user's
CV was too easy (optimistic); negative means it was pessimistic (e.g.
over-blocked) – both are reported faithfully.

## Usage

``` r
estimate_optimism(
  data,
  split,
  response,
  predict_fun = NULL,
  coords = NULL,
  metric = c("rmse", "mae", "brier", "logloss", "poisson"),
  control = c("block", "buffer"),
  dependence = NULL,
  k_control = NULL,
  buffer = NULL
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

- predict_fun:

  A function `(train_df, test_df)` returning numeric predictions for the
  test rows. Defaults to inverse-distance weighting on the coordinates.
  Custom learners (e.g. wrapping `lm`/`ranger`) let optimism be
  model-conditional.

- coords:

  For non-`sf` input, the coordinate column names/indices.

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

- dependence:

  Optional
  [`estimate_dependence()`](https://olatunjijohnson.github.io/spLeakage/reference/estimate_dependence.md)
  object (used to set the buffer for `control = "buffer"`).

- k_control:

  Number of control folds (defaults to the number of user folds).

- buffer:

  Buffer distance for `control = "buffer"` (defaults to the practical
  range from `dependence`).

## Value

An object of class `optimism_estimate`.

## Examples

``` r
set.seed(1)
d <- data.frame(x = runif(120), y = runif(120), z = rnorm(120))
estimate_optimism(d, split = sample(rep_len(1:10, 120)), response = "z",
                  coords = c("x", "y"))
#> <optimism_estimate>
#>   metric / control  : RMSE / block
#>   user CV error     : 1.028
#>   controlled error  : 1.039
#>   optimism          : +0.01121  (+1.1% of controlled)  [OPTIMISTIC (reported accuracy inflated)]
```
