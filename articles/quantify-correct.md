# Quantify, correct, and design

## From “is it leaking?” to “what do I do about it?”

Detecting leakage is only the first step. Once you know a split leaks,
four practical questions follow, and `spLeakage` answers each:

1.  **How inflated is my reported accuracy?** →
    [`estimate_optimism()`](https://olatunjijohnson.github.io/spLeakage/reference/estimate_optimism.md)
    (refit) and
    [`predict_optimism()`](https://olatunjijohnson.github.io/spLeakage/reference/predict_optimism.md)
    (fast emulator).
2.  **What number should I report instead?** →
    [`deleak_estimate()`](https://olatunjijohnson.github.io/spLeakage/reference/deleak_estimate.md)
    — a corrected accuracy with a confidence interval, even for a
    *published* result.
3.  **Which validation scheme should I have used?** →
    [`rank_cv_schemes()`](https://olatunjijohnson.github.io/spLeakage/reference/rank_cv_schemes.md).
4.  **Where should I collect independent validation data?** →
    [`design_validation()`](https://olatunjijohnson.github.io/spLeakage/reference/design_validation.md).

We reuse the clustered example throughout.

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

## 1. Quantify the optimism (the refit route)

[`estimate_optimism()`](https://olatunjijohnson.github.io/spLeakage/reference/estimate_optimism.md)
re-evaluates the model under a leakage-controlled, design-matched scheme
and reports the gap to your reported (random-CV) error.

``` r

estimate_optimism(d, folds, response = "z", coords = c("x", "y"))
#> <optimism_estimate>
#>   metric / control  : RMSE / block
#>   user CV error     : 0.2119
#>   controlled error  : 0.7568
#>   optimism          : +0.5449  (+72.0% of controlled)  [OPTIMISTIC (reported accuracy inflated)]
```

The `optimism` line is the headline: it gives the absolute gap *and* the
percentage inflation of your accuracy. A large positive percentage means
the random split made the model look that much better than the
design-matched scheme says it is. (For a true probability sample this
number can be negative — pessimism — and is reported as such.)

### Proper scoring rules for non-Gaussian responses

For prevalence (binary) or count responses, RMSE is not the right
metric. Pass a proper scoring rule so the optimism is meaningful —
`brier`/`logloss` for binary, `poisson` for counts:

``` r

db <- d; db$z <- rbinom(400, 1, plogis(2 * sin(2 * pi * d$x)))
estimate_optimism(db, folds, response = "z", coords = c("x", "y"), metric = "brier")
#> <optimism_estimate>
#>   metric / control  : BRIER / block
#>   user CV error     : 0.1798
#>   controlled error  : 0.1901
#>   optimism          : +0.0103  (+5.4% of controlled)  [OPTIMISTIC (reported accuracy inflated)]
```

## 2. Quantify the optimism (the fast emulator route)

The refit route needs to re-run the model.
[`predict_optimism()`](https://olatunjijohnson.github.io/spLeakage/reference/predict_optimism.md)
instead uses an emulator **pre-calibrated on a large simulation study**
and shipped with the package, returning an estimate (with a 90%
interval) in milliseconds and *without* refitting.

``` r

predict_optimism(d, folds, tgt, response = "z", coords = c("x", "y"), model = "idw")
#> Warning: Query is outside the emulator's area of applicability; refusing a
#> point estimate. Use estimate_optimism() instead.
#> <optimism_prediction> (emulator)
#>   model / response : idw / gaussian  [match: model+response]
#>   optimism         : NA  [outside AOA: DI 0.54 > 0.38]
```

Two things to note. First, when the query is in range the emulator
returns a point estimate with a genuine (conformal) predictive interval,
not a point guess. Second — and this is deliberate — the emulator
carries its **own** area-of-applicability guard. In fact this small,
tightly clustered example falls *outside* the envelope the emulator was
trained on, so it prints `optimism: NA [outside AOA: DI ... > ...]` and
**refuses to answer** rather than extrapolate. A leakage tool must not
commit the sin it polices. When the emulator abstains like this, fall
back to the refit route in section 1
([`estimate_optimism()`](https://olatunjijohnson.github.io/spLeakage/reference/estimate_optimism.md)),
which always applies.

## 3. Correct it: the de-leaked estimate

[`deleak_estimate()`](https://olatunjijohnson.github.io/spLeakage/reference/deleak_estimate.md)
is the number you should actually report: a bias-corrected (“de-leaked”)
accuracy with a bootstrap confidence interval. Pass your reported metric
value via `reported=` and it will also translate *that specific number*
into its corrected counterpart — so you can correct a published result
without the original model or refitting.

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

- `de-leaked estimate` (with CI) is the corrected accuracy.
- The `inflated by ...%` line restates the optimism as a multiplier on
  your error.
- The final line maps your supplied `reported` 0.30 to its de-leaked
  value.

For extrapolation deployments, pass `target=` so the correction uses a
**target-matched** control (the paper shows this recovers ~89% of the
true independent error, versus ~58% with a fixed block control).

## 4. Recommend and rank validation schemes

[`recommend_validation()`](https://olatunjijohnson.github.io/spLeakage/reference/recommend_validation.md)
gives design-aware guidance (recall: the design is *declared*, never
inferred):

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
```

Where
[`recommend_validation()`](https://olatunjijohnson.github.io/spLeakage/reference/recommend_validation.md)
gives expert guidance,
[`rank_cv_schemes()`](https://olatunjijohnson.github.io/spLeakage/reference/rank_cv_schemes.md)
makes it **data-driven**: it builds each candidate scheme on *your* data
and ranks them by how closely each one’s geometry matches your declared
deployment target (smallest `|SLI_rho|` wins).

``` r

rank_cv_schemes(d, tgt, response = "z", coords = c("x", "y"))
#> <scheme_ranking>  target = 'grid'
#>   ranked by deployment match (|SLI_rho| -> 0 is best):
#>     * block     SLI_rho +0.018  W 0.0254  [well matched]
#>       random    SLI_rho +0.093  W 0.0933  [OPTIMISTIC leakage]
#>       buffered  SLI_rho -0.197  W 0.222  [PESSIMISTIC (anti-leakage)]
#>   recommended: block
```

The `*` marks the recommended scheme. On clustered data you will
typically see random CV flagged optimistic, buffered CV *over-corrected
into pessimism*, and block CV closest to deployment — a nuance a static
lookup table cannot provide.

## 5. Design where to collect validation data

The honest way to know your true map accuracy is to measure it on
**independent** points. Given a budget, where should those points go?
[`design_validation()`](https://olatunjijohnson.github.io/spLeakage/reference/design_validation.md)
inverts the leakage logic: it stratifies space by distance-to-training
(on the deployment distribution), models the per-stratum error from the
spatial theory, and allocates the budget by Neyman allocation with
inclusion weights, so the resulting weighted-mean accuracy estimate is
both unbiased and minimum-variance.

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

The `SE ... lower than simple random validation` line quantifies the
saving: in the paper’s experiments the optimal design reaches the same
precision from a budget of 30 points that random placement needs far
more to match. The returned strata table (deployment weight × mean error
× points placed) shows where the budget went and the inclusion weights
to use when averaging the validation errors.

## Recap

`spLeakage` takes you from a single leakage score to a defensible,
corrected accuracy estimate and a plan for collecting the data that
would confirm it. Combined with the multi-channel detectors in
[`vignette("channels")`](https://olatunjijohnson.github.io/spLeakage/articles/channels.md)
and the conceptual framing in
[`vignette("spLeakage")`](https://olatunjijohnson.github.io/spLeakage/articles/spLeakage.md),
that is the full diagnostic workflow. \`\`\`
