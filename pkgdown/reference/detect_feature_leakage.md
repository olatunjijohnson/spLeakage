# Detect feature-space (covariate) leakage

Compares how close test points are to training points in standardised
covariate space against how close the deployment locations are to the
sample. A positive index means the test set is easier (closer in feature
space) than deployment – optimistic feature-space leakage / hidden
extrapolation at deployment.

## Usage

``` r
detect_feature_leakage(data, split, covariates, newdata = NULL, coords = NULL)
```

## Arguments

- data:

  An `sf` object or `data.frame` containing the covariate columns.

- split:

  A split specification (see
  [`detect_leakage()`](https://olatunjijohnson.github.io/spLeakage/reference/detect_leakage.md)).

- covariates:

  Character vector of covariate column names.

- newdata:

  Optional deployment data with the same covariate columns. Without it,
  only the test-set area-of-applicability is reported (no signed index).

- coords:

  Unused for the index (kept for API symmetry).

## Value

An object of class `feature_leakage`.
