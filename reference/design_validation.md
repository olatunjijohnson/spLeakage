# Design an optimal independent validation sample

Given training locations and a deployment target, choose where to place
a budget of independent validation points so the resulting estimate of
true map accuracy is unbiased (its distance-to-training distribution
matches deployment) and minimum-variance (distance-stratified Neyman
allocation using the GP error model).

## Usage

``` r
design_validation(
  data,
  target,
  budget,
  response = NULL,
  coords = NULL,
  candidates = NULL,
  allocation = c("optimal", "proportional", "random"),
  strata = 6L,
  dependence = NULL,
  seed = 1L
)
```

## Arguments

- data:

  Training data (`sf`, matrix, or `data.frame`).

- target:

  A [`prediction_target()`](prediction_target.md) (the deployment
  region).

- budget:

  Number of validation points to place.

- response:

  Name of the response column (for the dependence/error model).

- coords:

  For non-`sf` input, the coordinate column names/indices.

- candidates:

  Optional pool of candidate validation locations (default: the target
  locations).

- allocation:

  `"optimal"` (Neyman, default), `"proportional"`, or `"random"`.

- strata:

  Number of distance-to-training strata.

- dependence:

  Optional [`estimate_dependence()`](estimate_dependence.md) object.

- seed:

  RNG seed.

## Value

An object of class `validation_design`.

## Examples

``` r
set.seed(1)
d <- data.frame(x = runif(80), y = runif(80), z = rnorm(80))
grid <- as.matrix(expand.grid(x = seq(0, 1, .1), y = seq(0, 1, .1)))
design_validation(d, prediction_target(grid = grid, type = "grid"),
                  budget = 20, response = "z", coords = c("x", "y"))
#> <validation_design>  budget = 20, allocation = optimal, placed = 18
#>   estimate true accuracy as a weighted mean of validation errors (weights returned)
#>   SE of the accuracy estimate: 0% lower than simple random validation (same budget)
#>   distance strata (deployment weight x mean error x points placed):
#>  weight mean_error n
#>   0.174          1 3
#>   0.165          1 3
#>   0.165          1 3
#>   0.165          1 3
#>   0.157          1 3
#>   0.174          1 3
```
