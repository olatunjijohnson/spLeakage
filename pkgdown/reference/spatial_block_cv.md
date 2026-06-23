# Spatial block cross-validation folds

Assigns observations to `k` spatially contiguous folds by k-means
clustering of the coordinates. A leakage-controlled scheme for
comparison against a user's (often random) split.

## Usage

``` r
spatial_block_cv(data, k = 10L, coords = NULL)
```

## Arguments

- data:

  An `sf` object, numeric coordinate matrix, or `data.frame`.

- k:

  Number of spatial folds.

- coords:

  For non-`sf` input, the coordinate column names/indices.

## Value

An integer fold-id vector of length `nrow(data)`.
