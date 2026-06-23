# Audit a cross-validation workflow for spatial leakage

Runs
[`detect_leakage()`](https://olatunjijohnson.github.io/spLeakage/reference/detect_leakage.md)
and adds data-hygiene checks (duplicated coordinates, documented CRS) to
produce a scorecard.

## Usage

``` r
audit_workflow(
  data,
  split,
  target,
  response = NULL,
  dependence = NULL,
  coords = NULL,
  group = NULL
)
```

## Arguments

- data:

  An `sf` object, numeric coordinate matrix, or `data.frame`.

- split:

  A split specification: a `list(test=, train=)`, a list of test-index
  vectors (k folds), or a fold-id vector of length `nrow(data)`.

- target:

  A
  [`prediction_target()`](https://olatunjijohnson.github.io/spLeakage/reference/prediction_target.md)
  (or coordinates coercible to one).

- response:

  Name of the numeric response column, used to estimate spatial
  dependence (required unless `dependence` is supplied).

- dependence:

  A precomputed
  [`estimate_dependence()`](https://olatunjijohnson.github.io/spLeakage/reference/estimate_dependence.md)
  object (optional; overrides `response`). Carries any anisotropy used.

- coords:

  For non-`sf` input, the coordinate column names/indices.

- group:

  Optional grouping column for
  [`detect_group_leakage()`](https://olatunjijohnson.github.io/spLeakage/reference/detect_group_leakage.md)
  (defaults to grouping by coordinates).

## Value

An object of class `workflow_audit`.
