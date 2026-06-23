# Leakage channels beyond distance

## Why “distance” is not the whole story

The geographic Spatial Leakage Index (the [getting-started
vignette](https://olatunjijohnson.github.io/spLeakage/articles/spLeakage.md))
catches the most common leak: test points that are close to training
points *in space*. But a train/test split can leak information through
several other channels, each of which inflates your accuracy in a way
the geographic index alone will not see. This vignette walks through
each channel, shows the `spLeakage` diagnostic for it, and — where one
exists — the paired fix.

A small clustered dataset serves as a running example.

``` r

set.seed(1)
nc <- 8; ce <- cbind(runif(nc, .1, .9), runif(nc, .1, .9)); cl <- sample(nc, 400, TRUE)
xy <- ce[cl, ] + cbind(rnorm(400, 0, .04), rnorm(400, 0, .04))
d <- data.frame(x = pmin(pmax(xy[, 1], 0), 1), y = pmin(pmax(xy[, 2], 0), 1))
d$z <- sin(2 * pi * d$x) + cos(2 * pi * d$y) + rnorm(400, 0, .15)
grid <- as.matrix(expand.grid(x = seq(0, 1, .07), y = seq(0, 1, .07)))
tgt <- prediction_target(grid = grid, type = "grid")
```

## 1. Geographic (autocorrelation) leakage

This is the core channel, recapped here for completeness: the signed
Spatial Leakage Index from the distance between the CV and deployment
nearest-neighbour distributions. A positive `SLI_rho` means optimistic
leakage.

``` r

detect_leakage(d, split = sample(rep_len(1:10, 400)), target = tgt, response = "z",
               coords = c("x", "y"))
#> <leakage_diagnosis>
#>   target            : grid   |  n = 400, test = 400, folds = 10
#>   SLI_rho (signed)  : +0.092   [OPTIMISTIC leakage]
#>   SLI_d  (signed)   : +0.034   (A = +0.09289, phi = 2.72)
#>   retained corr.    : c_obs = 0.985 vs c_pred = 0.893
#>   W (NNDM) / delta  : 0.09289 / +1.00
```

See
[`vignette("spLeakage")`](https://olatunjijohnson.github.io/spLeakage/articles/spLeakage.md)
for a full reading of this output.

## 2. Grouped / duplicated-location leakage

A subtle but **exact** leak: the same site measured more than once
(repeat surveys, revisited plots, co-located records). If two records
from one location land in different folds, the model is tested on a
point it has, in effect, already seen.

We inject 60 duplicated locations and check a random split:

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
```

`test leaked via shared group` is the count (and fraction) of test
points whose location also appears in the training set. The fix,
[`group_kfold()`](https://olatunjijohnson.github.io/spLeakage/reference/group_kfold.md),
keeps all records from a location **together** in the same fold, driving
that fraction to zero:

``` r

gk <- group_kfold(dd, k = 10, coords = c("x", "y"))
detect_group_leakage(dd, gk, coords = c("x", "y"))$frac_leaked
#> [1] 0
```

This channel is what the malaria case study in the paper traced a
spurious geographic SLI back to: co-located repeat surveys, not genuine
autocorrelation.

## 3. Feature-space (covariate) leakage

A test point can be far in *geographic* space yet close in *covariate*
space — and then it is still easy to predict. Conversely, if deployment
requires extrapolating a covariate beyond its training range, your CV
(which never extrapolates) is optimistic. Supplying the deployment
covariates yields a signed feature-space index (an area-of-applicability
comparison, after Meyer & Pebesma 2021).

``` r

d$cov <- d$z + rnorm(400, 0, .3)
newdata <- data.frame(cov = runif(300, 1.2, 2.2))   # deployment extrapolates the covariate
detect_feature_leakage(d, sample(rep_len(1:5, 400)), covariates = "cov",
                       newdata = newdata)
#> <feature_leakage>
#>   covariates        : cov
#>   test in feature-AOA : 91% (interpolation in covariate space)
#>   feature SLI (signed): +3.816  [OPTIMISTIC feature leakage]
#>   mean NN reach (rel) : test 1.22 vs deployment 5.04
```

`test in feature-AOA` reports how much of the test set is interpolation
in covariate space; the signed `feature SLI` compares CV’s
covariate-space reach to deployment’s. A positive value warns that
deployment will go where CV never tested.

## 4. Covariate-extraction-overlap leakage

Common in remote sensing and SDMs: covariates are extracted from a
*window* (a focal buffer of radius `r`) around each point. If a test
window overlaps a training window, they share input pixels — leakage
through the feature construction itself, before any model is fit.

``` r

detect_extraction_leakage(d, sample(rep_len(1:5, 400)), radius = 0.08,
                          coords = c("x", "y"))
#> <extraction_leakage>  extraction radius = 0.08
#>   test windows overlapping a training window (< 2r): 400 / 400 (100.0%)
#>   test windows containing a training location (< r): 400 / 400 (100.0%)
#>   [!] fix: separate train/test by a buffer >= 0.16 (= 2 x radius)
```

The diagnostic counts test windows that overlap a training window
(closer than `2r`) or that contain a training location (closer than
`r`), and recommends a separating buffer of `2r`.

## 5. Temporal and joint spatiotemporal leakage

With a time dimension, a random split can train on the **future** to
predict the past — impossible at deployment. We build a small space–time
panel (each location revisited over 10 time steps):

``` r

st <- d[rep(1:40, each = 10), ]; st$t <- rep(1:10, 40)
detect_temporal_leakage(st, sample(rep_len(1:10, nrow(st))), time = "t")
#> <temporal_leakage>
#>   time column       : t   |  n test = 400
#>   lookahead leakage : 360 / 400 (90.0%) test points trained on the future
#>   mean time-to-nearest-training : 0
#>   [!] fix: temporal_kfold() (forward-chaining CV)
```

`lookahead leakage` is the fraction of test points whose nearest
training record lies in their future. The forward-chaining fix
([`temporal_kfold()`](https://olatunjijohnson.github.io/spLeakage/reference/temporal_kfold.md))
removes it:

``` r

fc <- temporal_kfold(st, k = 5, time = "t")
detect_temporal_leakage(st, fc, time = "t")$lookahead_frac
#> [1] 0
```

Crucially, leakage can hide in the *joint* space–time structure even
when the space-only and time-only views each look acceptable.
[`detect_st_leakage()`](https://olatunjijohnson.github.io/spLeakage/reference/detect_st_leakage.md)
measures all three in a single scaled space–time metric:

``` r

detect_st_leakage(st, sample(rep_len(1:10, nrow(st))), time = "t", coords = c("x", "y"))
#> <st_leakage>  (one-step-ahead forecast deployment)
#>   space-only SLI : +0.000  [matched]
#>   time-only  SLI : +0.370  [optimistic leakage]
#>   JOINT space-time SLI : +0.228  [optimistic leakage]
#>   median space-time NN: test 0.115 vs deployment 0.370 (in dependence units)
```

If the `JOINT space-time SLI` is larger than either marginal, a separate
spatial or temporal analysis would have missed the leak.

## 6. Trend strength (a leak the SLI cannot see)

The meta-audit in the paper turned up a second source of optimism that
**no** autocorrelation-based diagnostic (NNDM, the SLI, the de-leak) can
detect: a large-scale **trend** a trend-blind model cannot extrapolate.
Smooth fields like rainfall or elevation can be 40–54% optimistic yet
show almost no short-range leakage.
[`trend_strength()`](https://olatunjijohnson.github.io/spLeakage/reference/trend_strength.md)
flags this directly as the variance a quadratic coordinate surface
explains.

``` r

trendy <- data.frame(x = runif(200), y = runif(200))
trendy$z <- 3 * trendy$x - 2 * trendy$y + rnorm(200, 0, .2)
trend_strength(trendy, "z", coords = c("x", "y"))
#> [1] 0.9633468
```

A value near 1 means the response is dominated by a smooth coordinate
trend — treat extrapolation accuracy with suspicion regardless of what
the SLI says.

## Audit every channel at once

[`audit_workflow()`](https://olatunjijohnson.github.io/spLeakage/reference/audit_workflow.md)
runs the geographic, grouped, trend and CRS checks together and prints a
one-screen scorecard with `[ok]`/`[!]` flags — the quickest way to
triage a new analysis.

``` r

audit_workflow(dd, split = sample(rep_len(1:10, nrow(dd))), target = tgt,
               response = "z", coords = c("x", "y"))
#> <workflow_audit>
#>   leakage grade     : B  (SLI_rho = +0.095)
#>   [ok] autocorrelation leakage (random-split optimism)
#>   [!] trend strength: 0.95 (extrapolation-optimism risk a trend-blind model misses)
#>   [!] duplicated coordinates: 60
#>   [!] grouped/co-location leakage: 23.9% of test points
#>   [!] CRS documented: FALSE
```

To turn detected leakage into a number (“how inflated is my accuracy?”)
and to correct it, continue with
[`vignette("quantify-correct")`](https://olatunjijohnson.github.io/spLeakage/articles/quantify-correct.md).
\`\`\`
