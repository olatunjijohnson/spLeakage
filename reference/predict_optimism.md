# Predict optimism from the calibrated emulator (fast, no model refitting)

Estimates the optimism of a split using the pre-trained emulator and
cheap data/split/target features – the fast counterpart to
[`estimate_optimism()`](estimate_optimism.md). The estimate is
conditioned on the model class and response type. If the query falls
outside the emulator's area of applicability (in feature or category
space) the function refuses a point estimate (`NA`) rather than
extrapolating.

## Usage

``` r
predict_optimism(
  data,
  split,
  target,
  response = NULL,
  model = c("idw", "rf", "gam"),
  response_type = c("auto", "gaussian", "poisson", "binomial"),
  dependence = NULL,
  coords = NULL,
  emulator = NULL
)
```

## Arguments

- data:

  An `sf` object, numeric coordinate matrix, or `data.frame`.

- split:

  A split specification: a `list(test=, train=)`, a list of test-index
  vectors (k folds), or a fold-id vector of length `nrow(data)`.

- target:

  A [`prediction_target()`](prediction_target.md) (or coordinates
  coercible to one).

- response:

  Name of the numeric response column, used to estimate spatial
  dependence (required unless `dependence` is supplied).

- model:

  The learner whose optimism to estimate: `"idw"` (default), `"rf"`, or
  `"gam"`.

- response_type:

  `"auto"` (detect from the response, default), `"gaussian"`,
  `"poisson"`, or `"binomial"`.

- dependence:

  A precomputed [`estimate_dependence()`](estimate_dependence.md) object
  (optional; overrides `response`). Carries any anisotropy used.

- coords:

  For non-`sf` input, the coordinate column names/indices.

- emulator:

  An `optimism_emulator`; defaults to the one shipped with the package.

## Value

An object of class `optimism_prediction`.
