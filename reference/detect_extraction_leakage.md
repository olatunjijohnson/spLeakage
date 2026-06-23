# Detect covariate-extraction-overlap leakage

Flags test points whose covariate-extraction window overlaps a training
point's, given the extraction radius (focal-window radius, buffer
distance, or kernel bandwidth) used to build the spatial covariates.

## Usage

``` r
detect_extraction_leakage(data, split, radius, coords = NULL)
```

## Arguments

- data:

  An `sf` object, numeric matrix, or `data.frame`.

- split:

  A split specification (see
  [`detect_leakage()`](https://olatunjijohnson.github.io/spLeakage/reference/detect_leakage.md)).

- radius:

  The extraction-window radius / buffer / bandwidth, in coordinate units
  (geodesic metres for geographic CRS).

- coords:

  For non-`sf` input, the coordinate column names/indices.

## Value

An object of class `extraction_leakage`.

## Examples

``` r
set.seed(1)
d <- data.frame(x = runif(100), y = runif(100), z = rnorm(100))
# covariates were extracted with a 0.1-unit focal window:
detect_extraction_leakage(d, split = sample(rep_len(1:5, 100)),
                          radius = 0.1, coords = c("x", "y"))
#> <extraction_leakage>  extraction radius = 0.1
#>   test windows overlapping a training window (< 2r): 99 / 100 (99.0%)
#>   test windows containing a training location (< r): 91 / 100 (91.0%)
#>   [!] fix: separate train/test by a buffer >= 0.2 (= 2 x radius)
```
