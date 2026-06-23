# Detect temporal (lookahead) leakage in a split

Flags test observations whose fold's training set contains *later* time
points (the model is trained on the future to predict the past) and
summarises how close in time test points are to their training data.

## Usage

``` r
detect_temporal_leakage(data, split, time, coords = NULL)
```

## Arguments

- data:

  An `sf` object or `data.frame`.

- split:

  A split specification (see
  [`detect_leakage()`](https://olatunjijohnson.github.io/spLeakage/reference/detect_leakage.md)).

- time:

  Name of the time column (numeric or `Date`).

- coords:

  Unused (kept for API symmetry).

## Value

An object of class `temporal_leakage`.
