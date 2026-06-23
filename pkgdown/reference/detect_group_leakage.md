# Detect grouped / duplicated-location leakage in a split

Flags test observations that share a group (e.g. the same site, repeated
survey, household, or plot) with a training observation in the same fold
– exact leakage that the distance-based
[`detect_leakage()`](https://olatunjijohnson.github.io/spLeakage/reference/detect_leakage.md)
index registers only at distance zero.

## Usage

``` r
detect_group_leakage(data, split, group = NULL, coords = NULL, tol = 0)
```

## Arguments

- data:

  An `sf` object, numeric matrix, or `data.frame`.

- split:

  A split specification (see
  [`detect_leakage()`](https://olatunjijohnson.github.io/spLeakage/reference/detect_leakage.md)).

- group:

  Optional name of a grouping column. If omitted, groups are derived
  from the coordinates.

- coords:

  For non-`sf` input, the coordinate column names/indices.

- tol:

  When grouping by coordinates, the distance within which points are
  treated as the same location (`0` = exact duplicates).

## Value

An object of class `group_leakage`.

## Examples

``` r
# Two sites, each measured twice; a random split co-locates train and test.
d <- data.frame(x = c(1, 1, 2, 2), y = c(1, 1, 2, 2), z = rnorm(4))
detect_group_leakage(d, split = c(1, 2, 1, 2), coords = c("x", "y"))
#> <group_leakage>
#>   grouping        : coordinates (tol = 0)  |  n = 4, groups = 2, multi-member = 2
#>   test leaked via shared group : 4 / 4 (100.0%)
#>   groups split across folds    : 2
#>   [!] fix: group_kfold() keeps co-located/grouped records together
```
