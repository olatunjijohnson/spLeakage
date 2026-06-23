# Leakage channels beyond distance

Spatial autocorrelation is only one way a train/test split can leak.
`spLeakage` provides a diagnostic for each channel, and a paired fix
where one exists.

## 1. Geographic (autocorrelation) leakage

The core diagnostic: the signed Spatial Leakage Index.

``` r

set.seed(1)
nc <- 8; ce <- cbind(runif(nc, .1, .9), runif(nc, .1, .9)); cl <- sample(nc, 400, TRUE)
xy <- ce[cl, ] + cbind(rnorm(400, 0, .04), rnorm(400, 0, .04))
d <- data.frame(x = pmin(pmax(xy[, 1], 0), 1), y = pmin(pmax(xy[, 2], 0), 1))
d$z <- sin(2 * pi * d$x) + cos(2 * pi * d$y) + rnorm(400, 0, .15)
grid <- as.matrix(expand.grid(x = seq(0, 1, .07), y = seq(0, 1, .07)))
tgt <- prediction_target(grid = grid, type = "grid")

detect_leakage(d, split = sample(rep_len(1:10, 400)), target = tgt, response = "z",
               coords = c("x", "y"))
#> <leakage_diagnosis>
#>   target            : grid   |  n = 400, test = 400, folds = 10
#>   SLI_rho (signed)  : +0.092   [OPTIMISTIC leakage]
#>   SLI_d  (signed)   : +0.034   (A = +0.09289, phi = 2.72)
#>   retained corr.    : c_obs = 0.985 vs c_pred = 0.893
#>   W (NNDM) / delta  : 0.09289 / +1.00
```

## 2. Grouped / duplicated-location leakage

Repeated measurements at the same site, split across folds, leak
exactly.

``` r

dd <- rbind(d, d[sample(400, 60), ])            # 60 duplicated locations
gl <- detect_group_leakage(dd, split = sample(rep_len(1:10, nrow(dd))),
                           coords = c("x", "y"))
gl
#> <group_leakage>
#>   grouping        : coordinates (tol = 0)  |  n = 460, groups = 400, multi-member = 60
#>   test leaked via shared group : 112 / 460 (24.3%)
#>   groups split across folds    : 56
#>   [!] fix: group_kfold() keeps co-located/grouped records together
# fix: keep co-located records together
gk <- group_kfold(dd, k = 10, coords = c("x", "y"))
detect_group_leakage(dd, gk, coords = c("x", "y"))$frac_leaked
#> [1] 0
```

## 3. Feature-space (covariate) leakage

Test points close to training in *covariate* space are also easy. Supply
the deployment covariates to get a signed index.

``` r

d$cov <- d$z + rnorm(400, 0, .3)
newdata <- data.frame(cov = runif(300, 1.2, 2.2))   # deployment extrapolates the covariate
detect_feature_leakage(d, sample(rep_len(1:5, 400)), covariates = "cov", newdata = newdata)
#> <feature_leakage>
#>   covariates        : cov
#>   test in feature-AOA : 91% (interpolation in covariate space)
#>   feature SLI (signed): +3.816  [OPTIMISTIC feature leakage]
#>   mean NN reach (rel) : test 1.22 vs deployment 5.04
```

## 4. Covariate-extraction-overlap leakage

If covariates are built from a focal window of radius `r`, test windows
that overlap training windows leak. The fix is a buffer of `2r`.

``` r

detect_extraction_leakage(d, sample(rep_len(1:5, 400)), radius = 0.08, coords = c("x", "y"))
#> <extraction_leakage>  extraction radius = 0.08
#>   test windows overlapping a training window (< 2r): 400 / 400 (100.0%)
#>   test windows containing a training location (< r): 400 / 400 (100.0%)
#>   [!] fix: separate train/test by a buffer >= 0.16 (= 2 x radius)
```

## 5. Temporal and joint spatiotemporal leakage

``` r

st <- d[rep(1:40, each = 10), ]; st$t <- rep(1:10, 40)
# random CV trains on the future:
detect_temporal_leakage(st, sample(rep_len(1:10, nrow(st))), time = "t")
#> <temporal_leakage>
#>   time column       : t   |  n test = 400
#>   lookahead leakage : 360 / 400 (90.0%) test points trained on the future
#>   mean time-to-nearest-training : 0
#>   [!] fix: temporal_kfold() (forward-chaining CV)
# the joint space-time view catches leakage space-only analysis misses:
detect_st_leakage(st, sample(rep_len(1:10, nrow(st))), time = "t", coords = c("x", "y"))
#> <st_leakage>  (one-step-ahead forecast deployment)
#>   space-only SLI : +0.000  [matched]
#>   time-only  SLI : +0.370  [optimistic leakage]
#>   JOINT space-time SLI : +0.228  [optimistic leakage]
#>   median space-time NN: test 0.115 vs deployment 0.370 (in dependence units)
```

The fix for temporal leakage is forward-chaining CV:

``` r

fc <- temporal_kfold(st, k = 5, time = "t")
detect_temporal_leakage(st, fc, time = "t")$lookahead_frac
#> [1] 0
```

## 6. Trend strength (extrapolation risk)

A large-scale trend that a trend-blind model cannot extrapolate is a
leakage source the autocorrelation index does not see.

``` r

trendy <- data.frame(x = runif(200), y = runif(200))
trendy$z <- 3 * trendy$x - 2 * trendy$y + rnorm(200, 0, .2)
trend_strength(trendy, "z", coords = c("x", "y"))
#> [1] 0.9633468
```

All channels are summarised together by
[`audit_workflow()`](https://olatunjijohnson.github.io/spLeakage/reference/audit_workflow.md).
