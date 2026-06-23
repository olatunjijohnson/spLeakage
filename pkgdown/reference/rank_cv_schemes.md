# Rank candidate cross-validation schemes by deployment match

Builds each candidate CV scheme, measures its Spatial Leakage Index
against the declared prediction target, and ranks them: the scheme whose
`SLI_rho` is closest to zero best imitates deployment (neither
optimistic nor pessimistic).

## Usage

``` r
rank_cv_schemes(
  data,
  target,
  response,
  coords = NULL,
  k = 10L,
  schemes = c("random", "block", "buffered"),
  dependence = NULL,
  seed = 1L
)
```

## Arguments

- data:

  An `sf` object, numeric matrix, or `data.frame`.

- target:

  A [`prediction_target()`](prediction_target.md) (or coordinates
  coercible to one).

- response:

  Name of the numeric response column (for the dependence model).

- coords:

  For non-`sf` input, the coordinate column names/indices.

- k:

  Number of folds for each candidate scheme.

- schemes:

  Candidate schemes: any of `"random"`, `"block"`, `"buffered"`.

- dependence:

  Optional precomputed [`estimate_dependence()`](estimate_dependence.md)
  object.

- seed:

  RNG seed for the random scheme.

## Value

An object of class `scheme_ranking`.

## Examples

``` r
set.seed(1)
d <- data.frame(x = runif(150), y = runif(150), z = rnorm(150))
grid <- as.matrix(expand.grid(x = seq(0, 1, .2), y = seq(0, 1, .2)))
rank_cv_schemes(d, prediction_target(grid = grid, type = "grid"),
                response = "z", coords = c("x", "y"))
#> <scheme_ranking>  target = 'grid'
#>   ranked by deployment match (|SLI_rho| -> 0 is best):
#>     * random    SLI_rho +0.001  W 0.0124  [well matched]
#>       block     SLI_rho -0.009  W 0.0735  [well matched]
#>       buffered  SLI_rho -0.031  W 0.359  [PESSIMISTIC (anti-leakage)]
#>   recommended: random
```
