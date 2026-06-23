# Declare the prediction target (deployment geometry)

The Spatial Leakage Index is defined relative to where the model will
actually be used. This constructor declares that deployment geometry.

## Usage

``` r
prediction_target(
  data = NULL,
  grid = NULL,
  newdata = NULL,
  type = c("grid", "newdata", "interpolation"),
  coords = NULL,
  n = 5000L
)
```

## Arguments

- data:

  Sample data (used to derive an interpolation grid, and for CRS).

- grid, newdata:

  Explicit prediction locations (`sf`, matrix, or `data.frame`) for
  `type = "grid"` / `"newdata"`.

- type:

  One of `"grid"` (wall-to-wall mapping), `"newdata"` (supplied
  locations), or `"interpolation"` (unsampled locations within the
  sampled domain, generated from `data`).

- coords:

  For non-`sf` input, the coordinate column names/indices.

- n:

  Number of points to generate for `type = "interpolation"`.

## Value

An object of class `prediction_target` holding the prediction
coordinates and geometry info.
