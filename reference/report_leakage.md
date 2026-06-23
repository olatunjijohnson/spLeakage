# Assemble a leakage report (scorecard)

Combines a leakage diagnosis with optional optimism and validation
recommendation into a single printable scorecard for journal submission.

## Usage

``` r
report_leakage(
  diagnosis,
  optimism = NULL,
  recommendation = NULL,
  deleak = NULL
)
```

## Arguments

- diagnosis:

  A
  [`detect_leakage()`](https://olatunjijohnson.github.io/spLeakage/reference/detect_leakage.md)
  result.

- optimism:

  An optional
  [`estimate_optimism()`](https://olatunjijohnson.github.io/spLeakage/reference/estimate_optimism.md)
  result.

- recommendation:

  An optional
  [`recommend_validation()`](https://olatunjijohnson.github.io/spLeakage/reference/recommend_validation.md)
  result.

- deleak:

  An optional
  [`deleak_estimate()`](https://olatunjijohnson.github.io/spLeakage/reference/deleak_estimate.md)
  result (the corrected accuracy).

## Value

An object of class `leakage_report`.
