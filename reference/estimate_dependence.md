# Estimate spatial dependence (variogram, range, correlation, anisotropy)

Fits an empirical semivariogram and an exponential model, returning the
correlation function `rho(h)`, the variogram range, the practical range
`phi`, optional geometric anisotropy, and the parameter covariance used
to propagate uncertainty into the Spatial Leakage Index.

## Usage

``` r
estimate_dependence(
  data,
  response,
  coords = NULL,
  n_bins = 15L,
  cutoff = NULL,
  anisotropy = NULL
)
```

## Arguments

- data:

  An `sf` object, numeric coordinate matrix, or `data.frame`.

- response:

  Name of the (numeric) response column used for the variogram.

- coords:

  For non-`sf` input, column names/indices of the coordinates.

- n_bins:

  Number of variogram bins.

- cutoff:

  Maximum lag distance; defaults to half the maximum pairwise distance.

- anisotropy:

  `NULL` (isotropic, default), `"auto"` (estimate geometric anisotropy
  from directional variograms; projected data only), or a
  `list(angle =, ratio =)` giving the major-axis angle (radians) and the
  minor/major range ratio in `(0, 1]`.

## Value

An object of class `sp_dependence` with elements `rho` (a function),
`range`, `practical_range` (`phi`), `psill`, `nugget`, `signal_prop`,
`anisotropy` (or `NULL`), `coef`, `vcov`, and `variogram`.
