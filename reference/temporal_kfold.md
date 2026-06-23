# Forward-chaining temporal cross-validation (the fix for lookahead leakage)

Splits the data into `k` time-ordered blocks and builds folds where each
block is tested using only earlier blocks for training – no future
information leaks into training. Returns a list of folds usable directly
by [`detect_leakage()`](detect_leakage.md) etc.

## Usage

``` r
temporal_kfold(data, k = 5L, time, coords = NULL)
```

## Arguments

- data:

  An `sf` object or `data.frame`.

- k:

  Number of time blocks.

- time:

  Name of the time column (numeric or `Date`).

- coords:

  Unused (kept for API symmetry).

## Value

A list of folds, each `list(test=, train=)`.
