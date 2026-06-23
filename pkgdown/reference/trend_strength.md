# Trend strength: large-scale structure that resists extrapolation

The proportion of response variance explained by a low-order polynomial
of the coordinates – a measure of large-scale spatial trend. High trend
strength means a trend-blind model will extrapolate poorly, a leakage
channel the autocorrelation [`detect_leakage()`](detect_leakage.md)
index does not capture.

## Usage

``` r
trend_strength(data, response, coords = NULL, degree = 2L)
```

## Arguments

- data:

  An `sf` object, numeric matrix, or `data.frame`.

- response:

  Name of the numeric response column.

- coords:

  For non-`sf` input, the coordinate column names/indices.

- degree:

  Polynomial degree of the coordinate trend surface (default 2).

## Value

A single number in `[0, 1]`: the trend `R^2`.

## Examples

``` r
set.seed(1)
d <- data.frame(x = runif(100), y = runif(100))
d$z <- 3 * d$x + rnorm(100, 0, 0.2)        # strong x-trend
trend_strength(d, "z", coords = c("x", "y"))
#> [1] 0.9466331
```
