# Detect and quantify spatial leakage in a train/test split

Audits an existing split or cross-validation scheme: computes the
Spatial Leakage Index (signed) in both a model-free distance form
(`SLI_d`) and a covariance-aware dependence form (`SLI_rho`), with a
per-point leakage decomposition. Positive values indicate optimistic
leakage (CV easier than deployment); negative values indicate pessimism.
See `docs/METHOD-SLI.md`.

## Usage

``` r
detect_leakage(
  data,
  split,
  target,
  response = NULL,
  dependence = NULL,
  coords = NULL,
  n_boot = 0L,
  k = 1L
)
```

## Arguments

- data:

  An `sf` object, numeric coordinate matrix, or `data.frame`.

- split:

  A split specification: a `list(test=, train=)`, a list of test-index
  vectors (k folds), or a fold-id vector of length `nrow(data)`.

- target:

  A [`prediction_target()`](prediction_target.md) (or coordinates
  coercible to one).

- response:

  Name of the numeric response column, used to estimate spatial
  dependence (required unless `dependence` is supplied).

- dependence:

  A precomputed [`estimate_dependence()`](estimate_dependence.md) object
  (optional; overrides `response`). Carries any anisotropy used.

- coords:

  For non-`sf` input, the coordinate column names/indices.

- n_boot:

  Number of Monte-Carlo draws from the variogram-fit uncertainty used to
  put a confidence interval on the SLI. `0` (default) skips it. The
  split geometry is exact, so only `rho`/`phi` uncertainty is
  propagated.

- k:

  Number of nearest training neighbours for the dependence-form
  `SLI_rho` (density-aware, noisy-OR retained correlation). `1`
  (default, recommended) is the single-nearest-neighbour index, which is
  near-sufficient for predicting optimism (cor ~ 0.98 with exact GP
  optimism). `k > 1` improves the magnitude match but reduces the
  correlation (it over-counts correlated neighbours), so it is offered
  as an option rather than the default.

## Value

An object of class `leakage_diagnosis`.

## Examples

``` r
set.seed(1)
d <- data.frame(x = runif(80), y = runif(80), z = rnorm(80))
folds <- sample(rep_len(1:5, 80))
grid <- as.matrix(expand.grid(x = seq(0, 1, 0.2), y = seq(0, 1, 0.2)))
tgt <- prediction_target(grid = grid, type = "grid")
detect_leakage(d, folds, tgt, response = "z", coords = c("x", "y"))
#> <leakage_diagnosis>
#>   target            : grid   |  n = 80, test = 80, folds = 5
#>   SLI_rho (signed)  : +0.000   [well matched]
#>   SLI_d  (signed)   : +0.026   (A = +0.005592, phi = 0.2176)
#>   retained corr.    : c_obs = 0.000 vs c_pred = 0.000
#>   W (NNDM) / delta  : 0.007047 / +0.79
```
