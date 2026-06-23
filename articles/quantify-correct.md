# Quantify, correct, and design

Detecting leakage is the first step. `spLeakage` also quantifies the
*consequence* (how inflated your accuracy is), corrects it, recommends a
better scheme, and designs where to collect independent validation data.

``` r

set.seed(1)
nc <- 8; ce <- cbind(runif(nc, .1, .9), runif(nc, .1, .9)); cl <- sample(nc, 400, TRUE)
xy <- ce[cl, ] + cbind(rnorm(400, 0, .04), rnorm(400, 0, .04))
d <- data.frame(x = pmin(pmax(xy[, 1], 0), 1), y = pmin(pmax(xy[, 2], 0), 1))
d$z <- sin(2 * pi * d$x) + cos(2 * pi * d$y) + rnorm(400, 0, .15)
folds <- sample(rep_len(1:10, 400))
grid <- as.matrix(expand.grid(x = seq(0, 1, .07), y = seq(0, 1, .07)))
tgt <- prediction_target(grid = grid, type = "grid")
```

## Quantify the optimism

How much is the reported accuracy inflated, relative to a
leakage-controlled scheme?

``` r

estimate_optimism(d, folds, response = "z", coords = c("x", "y"))
#> <optimism_estimate>
#>   metric / control  : RMSE / block
#>   user CV error     : 0.2119
#>   controlled error  : 0.7568
#>   optimism          : +0.5449  (+72.0% of controlled)  [OPTIMISTIC (reported accuracy inflated)]
```

For binary (disease-prevalence) or count responses, use a proper scoring
rule:

``` r

db <- d; db$z <- rbinom(400, 1, plogis(2 * sin(2 * pi * d$x)))
estimate_optimism(db, folds, response = "z", coords = c("x", "y"), metric = "brier")
#> <optimism_estimate>
#>   metric / control  : BRIER / block
#>   user CV error     : 0.1798
#>   controlled error  : 0.1901
#>   optimism          : +0.0103  (+5.4% of controlled)  [OPTIMISTIC (reported accuracy inflated)]
```

## Correct it: the de-leaked estimate

The de-leaked estimate is the accuracy you should report, with a
confidence interval. It can even correct a *published* number without
refitting.

``` r

deleak_estimate(d, folds, response = "z", coords = c("x", "y"),
                reported = 0.30, n_boot = 200)
#> <deleak_estimate>
#>   metric / control     : RMSE / block
#>   reported (your CV)   : 0.2119
#>   de-leaked estimate   : 0.7497   (90% CI 0.4034, 1.033)
#>   => reported accuracy inflated by +72% (error x3.54)
#>   your reported 0.3 -> de-leaked 1.062 (90% CI 0.5523, 1.472)
```

## Fast emulator route

[`predict_optimism()`](https://olatunjijohnson.github.io/spLeakage/reference/predict_optimism.md)
uses a calibrated emulator (no refitting); it refuses queries outside
its area of applicability.

``` r

predict_optimism(d, folds, tgt, response = "z", coords = c("x", "y"), model = "idw")
#> Warning: Query is outside the emulator's area of applicability; refusing a
#> point estimate. Use estimate_optimism() instead.
#> <optimism_prediction> (emulator)
#>   model / response : idw / gaussian  [match: model+response]
#>   optimism         : NA  [outside AOA: DI 0.54 > 0.38]
```

## Recommend a validation scheme

[`recommend_validation()`](https://olatunjijohnson.github.io/spLeakage/reference/recommend_validation.md)
gives design-aware guidance;
[`rank_cv_schemes()`](https://olatunjijohnson.github.io/spLeakage/reference/rank_cv_schemes.md)
ranks candidate schemes by how well they match deployment.

``` r

recommend_validation(estimand = "prediction", design = "clustered", target = "grid")
#> <validation_recommendation>
#>   estimand / design / target : prediction / clustered / grid
#>   spatial CV appropriate     : YES
#>   recommended:
#>     - NNDM / kNNDM CV
#>     - Spatial block CV
#>   avoid: Random k-fold CV (optimistic under spatial autocorrelation)
#>   rationale: Conditional predictive skill from a clustered sample for a 'grid' target: match the CV geometry to deployment (NNDM/kNNDM/buffered) so test points are as far from training as prediction points are from the sample.
rank_cv_schemes(d, tgt, response = "z", coords = c("x", "y"))
#> <scheme_ranking>  target = 'grid'
#>   ranked by deployment match (|SLI_rho| -> 0 is best):
#>     * block     SLI_rho +0.018  W 0.0254  [well matched]
#>       random    SLI_rho +0.093  W 0.0933  [OPTIMISTIC leakage]
#>       buffered  SLI_rho -0.197  W 0.222  [PESSIMISTIC (anti-leakage)]
#>   recommended: block
```

## Design where to validate

Given a budget of independent validation points, where should they go to
estimate true accuracy most precisely?

``` r

design_validation(d, tgt, budget = 25, response = "z", coords = c("x", "y"))
#> <validation_design>  budget = 25, allocation = optimal, placed = 24
#>   estimate true accuracy as a weighted mean of validation errors (weights returned)
#>   SE of the accuracy estimate: 25% lower than simple random validation (same budget)
#>   distance strata (deployment weight x mean error x points placed):
#>  weight mean_error n
#>   0.169      0.025 1
#>   0.164      0.069 1
#>   0.169      0.134 3
#>   0.164      0.213 4
#>   0.164      0.306 6
#>   0.169      0.427 9
```
