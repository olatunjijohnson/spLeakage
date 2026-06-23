# spLeakage

> Detect and quantify **spatial information leakage** in predictive
> modelling.

`spLeakage` is a *diagnostic* R package. Given spatial data and a
train/test split (or an existing modelling workflow), it:

- **detects and scores** how much a validation scheme leaks information
  through spatial autocorrelation (the **Spatial Leakage Index**),
- **estimates the optimism** that leakage induces in reported accuracy
  (*“your RMSE is likely ~14% optimistic”*), and
- **recommends an appropriate validation strategy** for your sampling
  design and prediction target — including when spatial cross-validation
  is *not* the right choice.

It **audits rather than generates** folds, interoperating with
established spatial cross-validation packages (CAST, blockCV,
spatialsample, mlr3spatiotempcv) rather than competing with them.

## Status

In development (phase P2 — diagnostic, optimism, and recommendation
engine implemented; `R CMD check` clean). See
[`docs/VISION.md`](docs/VISION.md) for the methodology, full feature
catalogue, and roadmap, [`docs/METHOD-SLI.md`](docs/METHOD-SLI.md) and
[`docs/METHOD-EMULATOR.md`](docs/METHOD-EMULATOR.md) for the locked
method specs, and [`docs/PAPER.md`](docs/PAPER.md) for the paper plan.

## Quick example

``` r

library(spLeakage)

# Declare where the model will actually be used (wall-to-wall map):
tgt <- prediction_target(grid = pred_grid, type = "grid")

# 1. Diagnose an existing split (signed leakage index + uncertainty):
lk  <- detect_leakage(data, split = my_folds, target = tgt,
                      response = "y", n_boot = 300)
lk            # SLI_rho, 90% CI, NNDM W, verdict
plot(lk)      # observed-vs-deployment NN-distance ECDFs

# 2. Quantify the optimism it causes:
estimate_optimism(data, split = my_folds, response = "y", control = "block")

# 3. Get a design-aware recommendation (design is declared, not inferred):
recommend_validation(data, estimand = "prediction",
                     design = "clustered", target = "grid")
```

See the package vignette for a full worked example.

## Installation

Not yet on CRAN. Install the development version locally with
`devtools::install("spLeakage")`.

## License

MIT © Olatunji Johnson
