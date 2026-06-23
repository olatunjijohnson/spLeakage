# Recommend a validation strategy for a sampling design and prediction target

Operationalises the design x estimand x target framework: it asks for
the estimand and the sampling design (which cannot be inferred from
coordinates) and returns a ranked, conditional recommendation –
including when spatial CV is *not* appropriate. Any supplied data is
used only to compute a clustering risk flag, never to infer the design.

## Usage

``` r
recommend_validation(
  data = NULL,
  estimand = c("prediction", "population"),
  design = c("unknown", "probability", "clustered", "convenience"),
  target = c("grid", "interpolation", "newdata"),
  coords = NULL
)
```

## Arguments

- data:

  Optional `sf`/matrix/`data.frame`, used only for the clustering flag.

- estimand:

  `"prediction"` (conditional predictive skill at locations, default) or
  `"population"` (population-mean map accuracy over a region).

- design:

  `"unknown"` (default), `"probability"`, `"clustered"`, or
  `"convenience"`. This is an elicited fact about data collection.

- target:

  `"grid"` (wall-to-wall, default), `"interpolation"`, or `"newdata"`.

- coords:

  For non-`sf` input, the coordinate column names/indices.

## Value

An object of class `validation_recommendation`.

## Examples

``` r
# A clustered sample to be mapped wall-to-wall: spatial CV is appropriate.
recommend_validation(estimand = "prediction", design = "clustered", target = "grid")
#> <validation_recommendation>
#>   estimand / design / target : prediction / clustered / grid
#>   spatial CV appropriate     : YES
#>   recommended:
#>     - NNDM / kNNDM CV
#>     - Spatial block CV
#>   avoid: Random k-fold CV (optimistic under spatial autocorrelation)
#>   rationale: Conditional predictive skill from a clustered sample for a 'grid' target: match the CV geometry to deployment (NNDM/kNNDM/buffered) so test points are as far from training as prediction points are from the sample.
# A probability sample: random CV is correct, spatial CV is over-pessimistic.
recommend_validation(estimand = "prediction", design = "probability", target = "grid")
#> <validation_recommendation>
#>   estimand / design / target : prediction / probability / grid
#>   spatial CV appropriate     : NO
#>   recommended:
#>     - Random CV (unbiased for a probability sample)
#>   avoid: Forcing spatial CV (over-pessimistic here)
#>   rationale: For a probability sample, random CV gives unbiased predictive skill; spatial CV would be pessimistic.
```
