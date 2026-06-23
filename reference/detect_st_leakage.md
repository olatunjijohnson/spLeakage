# Detect joint spatiotemporal leakage

Measures leakage in a scaled space-time metric (each axis divided by its
dependence range) for a split, relative to a deployment target (by
default a one-step-ahead forecast). Reports the joint space-time index
together with the space-only and time-only indices, since the joint
measure can differ from either.

## Usage

``` r
detect_st_leakage(
  data,
  split,
  time,
  coords = NULL,
  sp_range = NULL,
  t_range = NULL,
  horizon = NULL
)
```

## Arguments

- data:

  An `sf` object, numeric matrix, or `data.frame`.

- split:

  A split specification (see [`detect_leakage()`](detect_leakage.md)).

- time:

  Name of the time column (numeric or `Date`).

- coords:

  For non-`sf` input, the coordinate column names/indices.

- sp_range, t_range:

  Spatial / temporal dependence ranges used to scale the two axes
  (defaults: 30% of the spatial extent / temporal span). Supply
  estimated ranges for a calibrated metric.

- horizon:

  Forecast horizon for the default deployment target (defaults to the
  median temporal spacing).

## Value

An object of class `st_leakage`.
