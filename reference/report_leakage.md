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

  A [`detect_leakage()`](detect_leakage.md) result.

- optimism:

  An optional [`estimate_optimism()`](estimate_optimism.md) result.

- recommendation:

  An optional [`recommend_validation()`](recommend_validation.md)
  result.

- deleak:

  An optional [`deleak_estimate()`](deleak_estimate.md) result (the
  corrected accuracy).

## Value

An object of class `leakage_report`.
