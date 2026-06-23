# Group-aware k-fold assignment (the fix for grouped leakage)

Assigns whole groups to folds so that members of a group (e.g.
co-located records) are never split across train and test. The remedy
for the leakage that
[`detect_group_leakage()`](https://olatunjijohnson.github.io/spLeakage/reference/detect_group_leakage.md)
flags.

## Usage

``` r
group_kfold(data, k = 10L, group = NULL, coords = NULL, tol = 0)
```

## Arguments

- data:

  An `sf` object, numeric matrix, or `data.frame`.

- k:

  Number of folds.

- group:

  Optional name of a grouping column. If omitted, groups are derived
  from the coordinates.

- coords:

  For non-`sf` input, the coordinate column names/indices.

- tol:

  When grouping by coordinates, the distance within which points are
  treated as the same location (`0` = exact duplicates).

## Value

An integer fold-id vector of length `nrow(data)`.
